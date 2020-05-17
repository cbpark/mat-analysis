{-# LANGUAGE BangPatterns      #-}
{-# LANGUAGE GADTs             #-}
{-# LANGUAGE MultiWayIf        #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Main where

import           MAT.Combinatorics                 (correctPairs)
import           MAT.Helper

import           HEP.Kinematics.Antler
-- hep-utilities
import           HEP.Data.LHCO

import           Codec.Compression.GZip            (decompress)
import qualified Data.ByteString.Char8             as C
import           Data.ByteString.Lazy.Char8        (ByteString)
import qualified Data.ByteString.Lazy.Char8        as BL
import           Data.Double.Conversion.ByteString

import           Pipes
import           Pipes.ByteString                  (fromLazy)
import qualified Pipes.Prelude                     as P
--
import           Control.Monad                     (forever, unless, when)
import           Data.List                         (sort, sortBy)
import           System.Environment                (getArgs)
import           System.Exit                       (die)
import           System.IO

-- import           Debug.Trace

main :: IO ()
main = do
    args <- getArgs
    let lenArg = length args
    unless (lenArg == 1 || lenArg == 2) $
        die "-- Usage: atLHCO <LHCO file gzipped> [output]"

    let lhcFile = head args
    putStrLn $ "-- The input LHCO file is " <> lhcFile <> "."
    events <- decompress <$> BL.readFile lhcFile

    let writeOutput h =
            runEffect $ getLHCOEvent fromLazy events
            -- >-> P.take 10
            >-> basicSelection
            >-> takeDLEvent
            >-> P.map (calcVar 0 80.379 173)
            >-> printVar h

    if lenArg == 1
        then writeOutput stdout
        else do let outfile = args !! 1
                withFile outfile WriteMode $ \h -> do
                    BL.hPutStrLn h header
                    putStrLn "-- Calculating the event variables ..."
                    writeOutput h

                putStrLn $ "-- ... Done!\n"
                    <> "-- " <> outfile <> " has been generated."

basicSelection :: MonadIO m => Pipe Event Event m ()
basicSelection = forever $ do
  Event {..} <- await

  when (nev `mod` 1000 == 0) $
      liftIO $ putStrLn ("-- processed " ++ show nev ++ " events.")

  yield $ Event { nev      = nev
                , photon   = filter basicSelection' photon
                , electron = filter basicSelection' electron
                , muon     = filter basicSelection' muon
                , tau      = filter basicSelection' tau
                , jet      = filter basicSelection' jet
                , bjet     = filter basicSelection' bjet
                , met      = met }

basicSelection' :: PhyObj a -> Bool
basicSelection' ObjPhoton {}                           = False
basicSelection' (ObjElectron (Track (eta', _, pt')) _) = abs eta' < 2.4 && pt' > 20
basicSelection' (ObjMuon (Track (eta', _, pt')) _ _)   = abs eta' < 2.4 && pt' > 20
basicSelection' ObjTau {}                              = False
basicSelection' (ObjJet (Track (eta', _, pt')) _ _ )   = abs eta' < 2.4 && pt' > 30
basicSelection' (ObjBjet (Track (eta', _, pt')) _ _ _) = abs eta' < 2.4 && pt' > 30
basicSelection' (ObjMet (_, pt'))                      = pt' > 40
basicSelection' ObjUnknown                             = False

data DLEvent = DLEvent { leptons :: [FourMomentum]
                       , jets    :: [FourMomentum]
                       , ptmiss  :: TransverseMomentum
                       } deriving Show

takeDLEvent :: Monad m => Pipe Event (Maybe DLEvent) m ()
takeDLEvent = forever $ do
    ev@Event {..} <- await
    let electrons = fourMomentum <$> electron
        muons     = fourMomentum <$> muon
        leptons'  = take 2 $ sortBy ptCompare (electrons <> muons)
        mll       = invariantMass leptons'
        jets      = fourMomentum <$> jet
        bjets     = fourMomentum <$> bjet
        jets'     = take 2 $ bjets <> jets
    yield $ if | length leptons' < 2 || length jets' < 2
                 || not (basicSelection' met) -> Nothing
               | null bjet                    -> Nothing
               | pt (head leptons') < 25      -> Nothing
               | mll > 76 && mll < 106        -> Nothing
               | mll < 20                     -> Nothing
               | otherwise -> Just $ DLEvent { leptons = leptons'
                                             , jets    = jets'
                                             , ptmiss  = missingET ev }

data Var = Var { _mAT0s  :: ![Double]
               , _mATs   :: ![Double]
               , _mMAOS  :: ![Double]
               , _mTtrue :: !Double
               , _mT2    :: !Double
               , _met    :: !Double
               } deriving Show

printVar :: MonadIO m => Handle -> Consumer (Maybe Var) m ()
printVar h = forever $ do
    vars <- await
    liftIO $ case vars of
                 Nothing       -> return ()
                 Just Var {..} -> C.hPutStrLn h $
                     C.unwords (map (\m -> toFixed 4 m <> "  ") _mAT0s)
                     <> C.unwords (map (\m -> toFixed 4 m <> "  ") _mATs)
                     <> C.unwords (map (\m -> toFixed 4 m <> "  ") _mMAOS)
                     <> toFixed 4 _mTtrue
                     <> "  " <> toFixed 4 _mT2
                     <> "  " <> toFixed 4 _met

calcVar :: Double  -- ^ invisible particle mass
        -> Double  -- ^ intermediate
        -> Double  -- ^ heavy
        -> Maybe DLEvent
        -> Maybe Var
calcVar m0 m1' m1 (Just ev@DLEvent {..}) = do
    pVs <- selectP m0 m1' m1 ev
    at <- mkAntler m0 m1 (visibles pVs)

    let pVtot = momentumSum pVs
        (qx, qy) = pxpy pVtot
        AT _ mAT0s = calcAT at qx qy 0 0

        mTtrue' = mTtrue at ptmiss
        met = norm ptmiss

    return $ case mATMAOS at qx qy ptmiss of
                 Nothing -> Var { _mAT0s  = mAT0s
                                , _mATs   = sort . concat . replicate 4 $ mAT0s
                                , _mMAOS  = [0, 0, 0, 0]
                                , _mTtrue = mTtrue'
                                , _mT2    = 0
                                , _met    = met }
                 Just (mATs, mMAOS, mT2) -> Var { _mAT0s  = mAT0s
                                                , _mATs   = mkLen 16 mATs
                                                , _mMAOS  = mkLen 4 mMAOS
                                                , _mTtrue = mTtrue'
                                                , _mT2    = mT2
                                                , _met    = met }
calcVar _ _ _ Nothing = Nothing

selectP :: Double -> Double -> Double -> DLEvent -> Maybe [FourMomentum]
selectP m0 m1' m1 DLEvent {..} =
    if m1' > m1
    then Nothing
    else do let blPairs0 = tupleToList <$> zip leptons jets
                !mMax = sqrt (m1 * m1 - m1' * m1')
            blPairs <- correctPairs blPairs0 ptmiss m1 m1' m0 mMax
            return $ momentumSum <$> blPairs
  where
    tupleToList :: (a, a) -> [a]
    tupleToList (a, b) = [a, b]

header :: ByteString
header = BL.pack $ "# " <>
         foldl1 (\v1 v2 -> v1 <> ", " <> v2)
         (zipWith (\n v -> "(" <> show n <> ") " <> v) ([1..] :: [Int])
             [ "mAT01", "mAT02", "mAT03", "mAT04"
             , "mAT11(maos)", "mAT12(maos)", "mAT13(maos)", "mAT14(maos)"
             , "mAT21(maos)", "mAT22(maos)", "mAT23(maos)", "mAT24(maos)"
             , "mAT31(maos)", "mAT32(maos)", "mAT33(maos)", "mAT34(maos)"
             , "mAT41(maos)", "mAT42(maos)", "mAT43(maos)", "mAT44(maos)"
             , "mMAOS1", "mMAOS2", "mMAOS3", "mMAOS4"
             , "mTtrue", "mT2", "MET" ])
