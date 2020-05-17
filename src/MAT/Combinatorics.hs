{-# LANGUAGE MultiWayIf #-}

module MAT.Combinatorics (correctPairs) where

import HEP.Kinematics
import HEP.Kinematics.Variable (mT2Symm, maosMomentaSymmetric)

import Data.Maybe              (isJust, isNothing)

-- import Debug.Trace

newtype Pair a = Pair { pair :: (a, a) } deriving Show

testVar1 :: (Pair FourMomentum, Pair FourMomentum) -> Double -> Maybe Double
testVar1 pairs m = if mPair > m then Nothing else Just mPair
  where mPair = testVar1' pairs

testVar1' :: (Pair FourMomentum, Pair FourMomentum) -> Double
testVar1' (pair1, pair2) = max (massP pair1) (massP pair2)
  where massP (Pair (p, q)) = invariantMass [p, q]

testVar2 :: (Pair FourMomentum, Pair FourMomentum)
         -> TransverseMomentum
         -> Double
         -> Double
         -> Maybe Double
testVar2 (Pair (p1, q1), Pair (p2, q2)) ptmiss mT mNu =
    let (v1, v2) = (p1 + q1, p2 + q2)
        mT2 = mT2Symm v1 v2 ptmiss mNu
    in if mT2 > mT then Nothing else Just mT2

testVar3 :: (Pair FourMomentum, Pair FourMomentum)
         -> TransverseMomentum
         -> Double
         -> Double
         -> Double
         -> Maybe Double
testVar3 (Pair (p1, q1), Pair (p2, q2)) ptmiss mT mW mNu = do
    let mT2WW = mT2Symm q1 q2 ptmiss mNu
        (chi1s, chi2s, _) =  maosMomentaSymmetric mT2WW q1 q2 ptmiss mW mNu
        (v1, v2) = (p1 + q1, p2 + q2)
        pT1s = (+ v1) <$> chi1s
        pT2s = (+ v2) <$> chi2s
        pTs = pT1s <> pT2s

    if null pTs
       then Nothing
       else do let mMAOS = map mass pTs
                   t3 = sum $ map (\m -> abs (m - mT)) mMAOS
               return t3

testVar4 :: (Pair FourMomentum, Pair FourMomentum)
         -> TransverseMomentum
         -> Double
         -> Double
         -> Double
         -> Maybe Double
testVar4 (Pair (p1, q1), Pair (p2, q2)) ptmiss mT mW mNu = do
    let (v1, v2) = (p1 + q1, p2 + q2)
        mT2tt = mT2Symm v1 v2 ptmiss mNu
        (chi1s, chi2s, _) =  maosMomentaSymmetric mT2tt v1 v2 ptmiss mT mNu
        pW1s = (+ q1) <$> chi1s
        pW2s = (+ q2) <$> chi2s
        pWs = pW1s <> pW2s

    if null pWs
        then Nothing
        else do let mMAOS = map mass pWs
                    t4 = sum $ map (\m -> abs (m - mW)) mMAOS
                return t4

testVars :: (Pair FourMomentum, Pair FourMomentum)
         -> TransverseMomentum
         -> Double
         -> Double
         -> Double
         -> Double
         -> (Maybe Double, Maybe Double, Maybe Double, Maybe Double)
testVars pairs ptmiss mY mW mX mVV =
    ( testVar1 pairs mVV
    , testVar2 pairs ptmiss mY mX
    , testVar3 pairs ptmiss mY mW mX
    , testVar4 pairs ptmiss mY mW mX )

pairsToList :: (Pair FourMomentum,  Pair FourMomentum) -> [[FourMomentum]]
pairsToList (Pair (p1, q1), Pair (p2, q2)) = [[p1, q1], [p2, q2]]

correctPairs :: [[FourMomentum]]
             -> TransverseMomentum
             -> Double  -- ^ m_Y (top)
             -> Double  -- ^ m_{intermediate} (W)
             -> Double  -- ^ m_X (neutrino)
             -> Double  -- ^ m_{vv}^max (m_{bl})
             -> Maybe [[FourMomentum]]
correctPairs [[p1, q1], [p2, q2]] ptmiss mT mW mNu mBLmax = do
    let pairs  = (Pair (p1, q1), Pair (p2, q2))
        pairs' = (Pair (p1, q2), Pair (p2, q1))

        -- (t1 , t2 , t3 , t4 ) = testVars pairs  ptmiss mT mW mNu mBLmax
        -- (t1', t2', t3', t4') = testVars pairs' ptmiss mT mW mNu mBLmax
        (t1 , t2 , _, _) = testVars pairs  ptmiss mT mW mNu mBLmax
        (t1', t2', _, _) = testVars pairs' ptmiss mT mW mNu mBLmax

        [listPairs, listPairs'] = pairsToList <$> [pairs, pairs']

    if | isNothing t1 && isJust    t1' -> return listPairs'
       | isJust    t1 && isNothing t1' -> return listPairs
       | isNothing t2 && isJust    t2' -> return listPairs'
       | isJust    t2 && isNothing t2' -> return listPairs
       {-
       | isNothing t4 && isJust    t4' -> return listPairs'
       | isJust    t4 && isNothing t4' -> return listPairs
       | length (catMaybes [t2, t3, t4, t2', t3', t4']) == 6 ->
             do t2Val <- t2
                t3Val <- t3
                t4Val <- t4
                t2Val' <- t2'
                t3Val' <- t3'
                t4Val' <- t4'

                let deltaT2, deltaT3, deltaT4 :: Int
                    deltaT2 = if t2Val - t2Val' > 0 then 1 else -1
                    deltaT3 = if t3Val - t3Val' > 0 then 1 else -1
                    deltaT4 = if t4Val - t4Val' > 0 then 1 else -1

                return $ if deltaT2 + deltaT3 + deltaT4 > 0
                         then listPairs'
                         else listPairs
       -}
       | otherwise -> return $ if testVar1' pairs > testVar1' pairs'
                               then listPairs'
                               else listPairs

correctPairs _ _ _ _ _ _ = Nothing
