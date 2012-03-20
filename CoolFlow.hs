import Data.Array.Unboxed
import Data.List
import Text.Printf

n = 100
timeStep = 0.1
numSteps = 1

indicies = [[0..n + 1], [1..n], [2..n-1]]
zeroMatrix = createMatrix (head indicies) (\ (_, _) -> 0.0)

type Matrix = UArray (Int, Int) Float
data Bound = XHard | YHard | Soft
data State = State !Matrix !Matrix !Matrix
applyNTimes :: Int -> (a -> a) -> a -> a
applyNTimes times f val = foldl' (\m _ -> f m) val [1..times]

showMatrix :: Matrix -> String
showMatrix m = foldl (++) "" [showRow i ++ "\n" | i <- head indicies] 
                where showRow i = foldl (++) "" [printf "%.3f " (m ! (i, j)) | j <- head indicies]
showState :: State -> String
showState (State _ _ d) = showMatrix d

createMatrix :: [Int] -> ((Int, Int) -> Float) -> Matrix
createMatrix ind f = array ((0, 0), (n + 1, n + 1)) [((i, j), f (i, j)) | j <- ind, i <- ind] // [(ix, 0) | ix <- other] 
			where other = [(i, j) | i <- head indicies, j <- head indicies] \\ [(i, j) | i <- ind, j <- ind]

add :: Matrix -> Matrix -> Matrix
add a b = createMatrix (head indicies) (\ix -> a ! ix + b ! ix)

addPecularSource :: State -> State
addPecularSource (State u v d) = State u v (d `add` createMatrix (head indicies) spike) 
                                        where spike ix = if ix == (n `div` 2, n `div` 2) then 0.1 else 0.0
advance :: State -> State
advance s = solveDensity . addPecularSource . solveVelocity $ s

main :: IO()
main = putStrLn . showState . applyNTimes numSteps advance $ State z z z where z = zeroMatrix

solveDensity :: State -> State
solveDensity (State u v d) = State u v (advect Soft (diffuse d diffusion) u v) where diffusion = 0.0
setCorners :: Matrix -> Matrix
setCorners a = a // [((n + 1, n + 1),	(a ! (n , n + 1) + a ! (n + 1, n)) / 2), 
	             ((n + 1, 0), 	(a ! (n , 0) + a ! (n + 1, 1)) / 2),
		     ((0, n + 1),	(a ! (0, n) + a ! (1, n + 1)) / 2),
		     ((0, 0),           (a ! (1 , 0) + a ! (0, 1)) / 2)]

diffuse :: Matrix -> Float -> Matrix
diffuse c0 diff = linearSolver c0 a (1 + 4 * a) where a = timeStep * diff * fromIntegral n * fromIntegral n 

setBoundary :: Bound -> Matrix -> Matrix
setBoundary b m = setCorners (setSides b m)

setSides :: Bound -> Matrix -> Matrix
setSides XHard m = setSides Soft m // ([((0, i), - m ! (0, i)) | i <- indicies !! 1] ++ [((n + 1, i), - m ! (n + 1, i)) | i <- indicies !! 1])
setSides YHard m = setSides Soft m // ([((i, 0), - m ! (i, 0)) | i <- indicies !! 1] ++ [((i, n + 1), - m ! (i, n + 1)) | i <- indicies !! 1])
setSides Soft m = m // ([((0, i), m ! (1, i)) | i <- indicies !! 1] ++ [((n + 1, i), m ! (n, i)) | i <- indicies !! 1] ++
			[((i, 0), m ! (i, 1)) | i <- indicies !! 1] ++ [((i, n + 1), m ! (i, n)) | i <- indicies !! 1])

linearSolver :: Matrix -> Float -> Float -> Matrix
linearSolver x0 a c = applyNTimes 20 iterFun zeroMatrix 
			where iterFun m = setBoundary Soft (createMatrix (indicies !! 1) (update m) ) 
			      update m (i, j) = (a * (m ! (i - 1, j) + m ! (i + 1, j) + m ! (i, j - 1) + m ! (i, j + 1)) + x0 ! (i, j)) / c

curl :: Matrix -> Matrix  -> Matrix
curl u v = createMatrix (indicies !! 1) (\(i, j) -> (u ! (i, j + 1) - u ! (i, j - 1) - v ! (i + 1, j) + v ! (i - 1, j)) / 2)

archimedes :: Matrix -> Matrix
archimedes d = createMatrix (head indicies) archCell 
		where (gravity, push) = (0.001, 0.005)
		      archCell (i, j) = gravity * d ! (i, j) + push * (d ! (i, j) - ambient)
		      ambient = sum [d ! (i, j) | i <- indicies !! 1, j <- indicies !! 1]

solveVelocity :: State -> State
solveVelocity (State u v d) = State u' v' d
				where (u', v') = project (advect YHard uProjected uProjected vProjected) (advect XHard vProjected uProjected vProjected)
				      (uProjected, vProjected) = project uDiffused vDiffused
				      uDiffused = diffuse uVorticityConfined viscosity
				      vDiffused = diffuse (vVorticityConfined `add` archimedes d) viscosity
				      (uVorticityConfined, vVorticityConfined) = vorticityConfimentForce u v

				      viscosity = 0.0
advectCell :: Matrix -> Matrix -> Matrix -> (Int, Int) -> Float
advectCell d0 du dv (i, j) = s0 * t0 * d0 ! (floor x, floor y) + s0 * t1 * d0 ! (floor x, floor y + 1) + 
			     s1 * t0 * d0 ! (floor x + 1, floor y) + s1 * t1 * d0 ! (floor x + 1, floor y + 1)
				where (s0, t0, s1, t1) = (1.0 - x + xi, 1.0 - y + yi, x - xi, y - yi)
				      (xi, yi) = (fromIntegral (floor x), fromIntegral (floor y))
			              x = sort [0.5 + fromIntegral n, 0.5, fromIntegral i - dt0 * du ! (i, j)] !! 1
		    	              y = sort [0.5 + fromIntegral n, 0.5, fromIntegral j - dt0 * dv ! (i, j)] !! 1
			              dt0 = timeStep * fromIntegral n

advect :: Bound -> Matrix -> Matrix -> Matrix -> Matrix
advect b d0 du dv = setBoundary b (createMatrix (indicies !! 1) (advectCell d0 du dv))

project :: Matrix -> Matrix -> (Matrix, Matrix)
project x y = (setBoundary YHard (createMatrix (indicies !! 1) projectCellX), setBoundary XHard (createMatrix (indicies !! 1) projectCellY)) 
                where projectCellY (i, j) = y ! (i, j) - 0.5 * fromIntegral n * (p ! (i, j + 1) - p ! (i, j - 1))
      		      projectCellX (i, j) = x ! (i, j) - 0.5 * fromIntegral n * (p ! (i + 1, j) - p ! (i - 1, j))
		      p = linearSolver (setBoundary Soft (createMatrix (indicies !! 1) computeDiv)) 1 4
		      computeDiv (i, j) = -0.5 * (x ! (i + 1, j) - x ! (i - 1, j) + y ! (i, j + 1) - y ! (i, j - 1)) / fromIntegral n

vorticityConfimentForce :: Matrix -> Matrix -> (Matrix, Matrix)
vorticityConfimentForce u v = (createMatrix (indicies !! 2) genX, createMatrix (indicies !! 2) genY) 
                                where genX ix = - c ! ix * dwdy ix / len ix
                                      genY ix = c ! ix * dwdx ix / len ix
                                      len  ix = sqrt(dwdx ix * dwdx ix + dwdy ix * dwdy ix) + eps
                                      dwdx (i, j) = (abs (c ! (i + 1, j)) - abs (c ! (i - 1, j))) / 2
                                      dwdy (i, j) = (abs (c ! (i, j + 1)) - abs (c ! (i, j - 1))) / 2
                                      c = curl u v
                                      eps = 1e-6

