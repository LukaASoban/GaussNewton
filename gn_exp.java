import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

import Jama.Matrix;

public class gn_exp {

	public static void main(String[] args) throws IOException {
		
		//Asks the user for the path of the .txt data file
		System.out.println("This is the Gauss-Newton solver for a Exponential"
				+ " function. Please type in the name of the text file containing"
				+ " a list of points. (ex. data.txt)");
		Scanner sc = new Scanner(System.in);
		String fileName = sc.nextLine();
		File file = new File(fileName);
		FileReader reader = new FileReader(file);
		BufferedReader bRead = new BufferedReader(reader);

		//initialize the data structures
		ArrayList<Double> xValues = new ArrayList<Double>();
		ArrayList<Double> yValues = new ArrayList<Double>();
		double[] beta = new double[3];
		int iter = 0;

		String line;
		while((line = bRead.readLine()) != null) {
			String[] temp = line.split(",");
			xValues.add(Double.parseDouble(temp[0]));
			yValues.add(Double.parseDouble(temp[1]));
		}

		Matrix b = new Matrix(beta, beta.length);

		//this is where the user inputs the initial guesses for beta (a,b,c)
		System.out.println("What is your initial guess for parameter a?");
		b.set(0, 0, sc.nextDouble());
		System.out.println("What is your initial guess for parameter b?");
		b.set(1, 0, sc.nextDouble());
		System.out.println("What is your initial guess for parameter c?");
		b.set(2, 0, sc.nextDouble());

		//asks the user for the number of iterations
		System.out.println("How many iterations do you want to run the"
				+ " Gauss-Newton algorithm?");
		iter = sc.nextInt();

		//the user may either choose to do QR Factorization by Householder reflections
		//or by Givens rotations
		System.out.println("Which QR Factorization method would you like to run?\n"
				+ "Householder or Givens. Type \"H\" or \"G\"");
		String method = sc.next();

		sc.close();
		bRead.close();

		//code splits depending on users choice of Householder or Givens
		int n = xValues.size();
		double[] res = new double[n];
		if (method.equalsIgnoreCase("H")) {

			for (int k = 0; k < iter; k++) {
				for (int i = 0; i < n; i++) {
					res[i] = (yValues.get(i) - (b.get(0, 0) * Math.pow(Math.E, (b.get(1, 0) * xValues.get(i))) + b.get(2, 0)));
				}
				Matrix residuals = new Matrix(res, res.length);

				//creates the Jacobian Matrix for QR factorization
				double[][] jacobian = new double[n][3];
				for (int i = 0; i < n; i++) {
					jacobian[i][0] = -(Math.pow(Math.E, (b.get(1, 0) * xValues.get(i))));
				}
				for (int i = 0; i < n; i++) {
					jacobian[i][1] = -(b.get(0, 0) * xValues.get(i) * Math.pow(Math.E, (b.get(1, 0) * xValues.get(i))));
				}
				for (int i = 0; i < n; i++) {
					jacobian[i][2] = -1;
				}
				Matrix j = new Matrix(jacobian);

				//sends the matrix to method for QR by Householder reflections
				qr_fact_househ qr = new qr_fact_househ(j);
				Matrix q = qr.getQ();
				Matrix r = qr.getR();
				
				System.out.println("Iteration " + (k + 1) + " gives Q and R as:");
				
				System.out.println("Q =");
				q.print(q.getColumnDimension(), 3);
				
				System.out.println("R =");
				r.print(r.getColumnDimension(), 3);
				
				//solves for x to find the new convergence
				b = b.minus(solveForX(q, r, residuals));

			}
			
			//print what it converges to
			System.out.print("Using the Householder method, this converges to:");
			b.print(b.getColumnDimension(), 3);

		} else if (method.equalsIgnoreCase("G")) {

			for (int k = 0; k < iter; k++) {
				for (int i = 0; i < n; i++) {
					res[i] = (yValues.get(i) - (b.get(0, 0) * Math.pow(Math.E, (b.get(1, 0) * xValues.get(i))) + b.get(2, 0)));
				}
				Matrix residuals = new Matrix(res, res.length);

				//creates the Jacobian Matrix for QR factorization
				double[][] jacobian = new double[n][3];
				for (int i = 0; i < n; i++) {
					jacobian[i][0] = -(Math.pow(Math.E, (b.get(1, 0) * xValues.get(i))));
				}
				for (int i = 0; i < n; i++) {
					jacobian[i][1] = -(b.get(0, 0) * xValues.get(i) * Math.pow(Math.E, (b.get(1, 0) * xValues.get(i))));
				}
				for (int i = 0; i < n; i++) {
					jacobian[i][2] = -1;
				}
				Matrix j = new Matrix(jacobian);

				//sends the matrix to method for QR by Givens rotations
				qr_fact_givens givens = new qr_fact_givens(j);
				
				givens.qrFactorization();
				
				System.out.println("Iteration " + (k + 1) + " gives Q and R as:");
				Matrix q = givens.getQ();
				System.out.println("Q =");
				q.print(q.getColumnDimension(), 3);
				
				Matrix r = givens.getR();
				System.out.println("R =");
				r.print(r.getColumnDimension(), 3);
				
				//solves for x to find the new convergence
				b = b.minus(solveForX(q, r, residuals));
			}
			
			//print what it converges to
			System.out.print("Using the Givens method, this converges to:");
			b.print(b.getColumnDimension(), 3);

		} else {
			System.out.print("Invalid input. Program will quit.");
		}

	}

	//Solves for x using back substitution
	private static Matrix solveForX(Matrix q, Matrix r, Matrix residuals) {

		//sets the answers (Qt * r) to the system of equations (b1, b2, b3)
		Matrix b = q.transpose().times(residuals);
		//Creates new matrix x for the solutions to the system of equations
		Matrix x = new Matrix(3,1);

		//The system of equations will always have the same form
		//This code solves for x
		x.set(2, 0, (b.get(2, 0) / r.get(2, 2)));
		x.set(1, 0, ((b.get(1, 0) - r.get(1, 2) * x.get(2, 0)) / r.get(1, 1)));
		double temp = r.get(0, 1) * x.get(1, 0);
		temp = b.get(0, 0) - temp + r.get(0, 2) * x.get(2, 0);
		x.set(0, 0, (temp / r.get(0, 0)));
		//returns the matrix x containing the solutions to the system of equations
		return x;
	}
}
