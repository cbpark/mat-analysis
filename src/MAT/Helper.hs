{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module MAT.Helper where

import           HEP.Kinematics.Antler

import           HEP.Kinematics
import           HEP.Kinematics.Vector.LorentzTVector (setXYT)

import           Data.ByteString                      (ByteString)
import qualified Data.ByteString.Char8                as C
import           Data.Double.Conversion.ByteString    (toExponential, toFixed)

data AT = AT { _deltaAT :: Double  -- ^ Delta_{AT}
             , _mATs    :: [Double]
             } deriving Show

showAT :: AT -> ByteString
showAT AT {..} = toExponential 6 _deltaAT
                 <> C.unwords (map (\m -> "  " <> toFixed 4 m) _mATs)

calcAT :: Antler
       -> Double  -- ^ - p_{x} component of the ISR
       -> Double  -- ^ - p_{y} component of the ISR
       -> Double  -- ^ a guess of the longitudinal momentum of the resonance
       -> Double  -- ^ the energy of the resonance
       -> AT
calcAT at qx qy qz e =
    let deltaATval = deltaAT at qx qy qz e
    in case mAT at qx qy qz of
           Nothing           -> AT deltaATval [0, 0, 0, 0]
           Just mATs         -> AT { _deltaAT  = deltaATval
                                   , _mATs = mkLen 4 mATs }

mkLen :: Num a => Int -> [a] -> [a]
mkLen n xs | length xs >= n  = take n xs
           | otherwise       = mkLen n $! xs <> [0]

mTtrue :: Antler -> TransverseMomentum -> Double
mTtrue Antler {..} ptmiss =
    let pChiT = setXYT (px ptmiss) (py ptmiss) (norm ptmiss)
    in transverseMass [_v1, _v2] pChiT
