import Jama.Matrix;

public class qr_fact_givens {

	//instance variables
	public Matrix A, R, Q;
	
	//the final values R and Q for QR Factorization
	public Matrix qr_R;
	public Matrix qr_Q;
	
	//constructor
	//takes in a matrix and sets A to the matrix
	public qr_fact_givens(Matrix matrix) {
		this.A = matrix;
	}
	
	public void qrFactorization() {
		
		//get the the number of rows and columns and set them to m and n respectively
		int m = A.getRowDimension();
		int n = A.getColumnDimension();
		
		
		//Setting matrix Q to be an identity matrix
		Matrix Q = new Matrix(m, m);
		for (int x = 0; x < m; x++) {
			Q.set(x, x, 1);
		}
		
		//Setting matrix R to be our matrix that is passed to the constructor (the Jacobian)
		Matrix R = A;
		
		double[] cs = new double[2];
		
		//will loop through all of the columns
		for (int j = 0; j <= n - 1; j++) {
			
			//re initialize variable i to be the number of rows - 1 so that it loops through the
			//rows each time
			int i = m - 1 ;
			while(i >= j + 1) {
				
				//create a new identity matrix
				Matrix G = new Matrix(m, m);
				for (int x = 0; x < m; x++) {
					G.set(x, x, 1);
				}
				
				//grabs the values to be the [x,y] vector then do Givens rotation on them
				//to find the cos and sin values
				cs = givensRotation(R.get(i - 1, j), R.get(i, j));
				
				//change the identity matrix to have cos sin
				//									-sin cos
				G.set(i, i, cs[0]);
				G.set(i - 1, i, -cs[1]);
				G.set(i - 1, i - 1, cs[0]);
				G.set(i, i - 1, cs[1]);
				
				
				//R = Gn * Gn-1 * ... * A
				R = G.transpose().times(R);
				//Q = Gn * Gn+1 * ...
				Q = Q.times(G);
				
				i--;
			}
			
		}
		
		//Make the R a 3x3 Matrix
		qr_R = R.getMatrix(0,2,0,2);
		
		//make the Q an Nx3 Matrix
		qr_Q = Q.getMatrix(0,Q.getRowDimension() - 1,0,2);
	}
	
	public double[] givensRotation(double x, double y) {
		
		//create a new 2 element array to hold cos and sin values
		// cs[0] is cos and cs[1] is sin
		double[] cs = new double[2];
		double temp;
		
		if(y == 0) {
			cs[0] = 1;
			cs[1] = 0;
		} else {
			if(Math.abs(y) > Math.abs(x)) {
				temp = x / y;
				cs[1] = (1 / Math.sqrt(1 + (temp * temp)));
				cs[0] = cs[1] * temp;
			} else {
				temp = y / x;
				cs[0] = (1 / Math.sqrt(1 + (temp * temp)));
				cs[1] = cs[0] * temp;
			}
		}
		return cs;
	}
	
	// Methods that get R and Q
	public Matrix getR() {
		return this.qr_R;
	}
	
	public Matrix getQ() {
		return this.qr_Q;
	}
	
}