{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Main where

import           MAT.Combinatorics                 (correctPairs)
import           MAT.Helper                        as MH

import           HEP.Kinematics.Antler
-- hep-utilities
import           HEP.Data.LHEF
import           HEP.Kinematics.Vector.TwoVector   (setXY)

import           Codec.Compression.GZip            (decompress)
import qualified Data.ByteString.Char8             as C
import           Data.ByteString.Lazy.Char8        (ByteString)
import qualified Data.ByteString.Lazy.Char8        as BL
import           Data.Double.Conversion.ByteString

import           Pipes
import           Pipes.ByteString                  (fromLazy)
import qualified Pipes.Prelude                     as P
--
import           Control.Monad                     (forever, unless)
import           Data.List                         (sort)
import           System.Environment                (getArgs)
import           System.Exit                       (die)
import           System.IO

-- import           Debug.Trace

main :: IO ()
main = do
    args <- getArgs
    let lenArg = length args
    unless (lenArg == 1 || lenArg == 2) $
        die "-- Usage: atLHEF2 <LHEF file gzipped> [output]"

    let lheFile = head args
    putStrLn $ "-- The input LHEF file is " <> lheFile <> "."
    events <- decompress <$> BL.readFile lheFile

    let writeOutput h =
            runEffect $ getLHEFEvent fromLazy events
            -- >-> P.map (calcVar 80.379 173.0 800)
            >-> P.map (calcVar 0 173.0 800)
            -- >-> P.take 10
            >-> printVar h

    if lenArg == 1
        then writeOutput stdout
        else do let outfile = args !! 1
                withFile outfile WriteMode $ \h -> do
                    BL.hPutStrLn h header
                    putStrLn $ "-- Calculating the event variables"
                        <> " (with combinatorial errors) ..."
                    writeOutput h

                putStrLn $ "-- ... Done!\n"
                    <> "-- " <> outfile <> " has been generated."

data Var = Var { _AT0         :: !AT
               , _mATs        :: ![Double]
               , _mMAOS       :: ![Double]
               , _mTtrue      :: !Double
               , _mT2         :: !Double
               , _correctPair :: !Int
               } deriving Show

printVar :: MonadIO m => Handle -> Consumer (Maybe Var) m ()
printVar h = forever $ do
    vars <- await
    liftIO $ case vars of
                 Nothing       -> return ()
                 Just Var {..} -> C.hPutStrLn h $
                     showAT _AT0
                     <> C.unwords (map (\m -> "  " <> toFixed 4 m) _mATs)
                     <> C.unwords (map (\m -> "  " <> toFixed 4 m) _mMAOS)
                     <> "  " <> toFixed 4 _mTtrue
                     <> "  " <> toFixed 4 _mT2
                     <> "  " <> C.pack (show _correctPair)

calcVar :: Double -> Double -> Double -> Event -> Maybe Var
calcVar m0 m1 m2 ps = do
    (pH, pVs, ptmiss, pairp) <- selectP ps
    at <- mkAntler m0 m1 (visibles pVs)
    let (qx, qy) = pxpy pH
        correctPair = if pairp then 1 else -1
        at0 = calcAT at qx qy 0 m2

    return $
        case mATMAOS at qx qy ptmiss of
            Nothing -> Var { _AT0         = at0
                           , _mATs        = sort . concat . replicate 4 $
                                            MH._mATs at0
                           , _mMAOS       = [0, 0, 0, 0]
                           , _mTtrue      = mTtrue at ptmiss
                           , _mT2         = 0
                           , _correctPair = correctPair }
            Just (mATs, mMAOS, mT2) ->
                Var { _AT0         = at0
                    , _mATs        = mkLen 16 mATs
                    , _mMAOS       = mkLen 4 mMAOS
                    , _mTtrue      = mTtrue at ptmiss
                    , _mT2         = mT2
                    , _correctPair = correctPair }

selectP :: Event
        -> Maybe (FourMomentum, [FourMomentum], TransverseMomentum, Bool)
selectP ev = do
    let topChild = particlesFrom topQuarks (eventEntry ev)
    if null topChild
        then Nothing
        else do let pH = momentumSum $ fourMomentum <$> concat topChild
                    pBLs = fmap fourMomentum <$>
                          (filter (not . isNeutrino) <$> topChild)
                    pNu = momentumSum (momentumSum . fmap fourMomentum <$>
                                       (filter isNeutrino <$> topChild))
                    ptmiss = setXY (px pNu) (py pNu)

                pBLpairs <- correctPairs pBLs ptmiss 173.0 80.379 0 153.173

                let pVs = momentumSum <$> pBLpairs
                return (pH, pVs, ptmiss, pBLs == pBLpairs)
  where
    topQuarks = ParticleType [6]
    isNeutrino = (`elem` neutrinos) . idOf

header :: ByteString
header = BL.pack $ "# " <>
         foldl1 (\v1 v2 -> v1 <> ", " <> v2)
         (zipWith (\n v -> "(" <> show n <> ") " <> v) ([1..] :: [Int])
             [ "deltaAT(0)", "mAT01", "mAT02", "mAT03", "mAT04"
             , "mAT11(maos)", "mAT12(maos)", "mAT13(maos)", "mAT14(maos)"
             , "mAT21(maos)", "mAT22(maos)", "mAT23(maos)", "mAT24(maos)"
             , "mAT31(maos)", "mAT32(maos)", "mAT33(maos)", "mAT34(maos)"
             , "mAT41(maos)", "mAT42(maos)", "mAT43(maos)", "mAT44(maos)"
             , "mMAOS1", "mMAOS2", "mMAOS3", "mMAOS4"
             , "mTtrue", "mT2", "correct pair" ])
