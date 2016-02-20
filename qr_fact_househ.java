import Jama.Matrix;

	public class qr_fact_househ {

		public Matrix matrix;
		public Matrix vect;

		public qr_fact_househ(Matrix matrix) {
			//creates new vector of zeros
			this.vect = new Matrix(1, matrix.getColumnDimension());
			this.matrix = matrix.copy();

			//loops through the columns
			for (int j = 0; j < this.matrix.getColumnDimension(); j++) {

				//loop to find the norm
				double norm = 0;
				for (int i = j; i < this.matrix.getRowDimension(); i++) {
					norm += Math.pow(this.matrix.get(i, j), 2);
				}
				norm = Math.sqrt(norm);

				//makes sure the sign of the norm is correct
				int sign = 1;
				if (this.matrix.get(j, j) < 0) {
					sign = -1;
				}

				double v1 = this.matrix.get(j, j) + sign * norm;
				double scalar = 1;

				for (int i = j + 1; i < this.matrix.getRowDimension(); i++) {
					this.matrix.set(i, j, (this.matrix.get(i, j) / v1));
					scalar += Math.pow(this.matrix.get(i, j), 2);
				}

				vect.set(0, j, (2 / scalar));

				//makes a new vector of zeros
				Matrix w = new Matrix(1, this.matrix.getColumnDimension());

				//grabs the significant vector for each new householder matrix
				for (int i = j; i < this.matrix.getColumnDimension(); i++) {
					w.set(0, i, this.matrix.get(j, i));
					for (int k = j + 1; k < this.matrix.getRowDimension(); k++) {
						if (i == j) {
							w.set(0, i, w.get(0, i) + this.matrix.get(k, j) * this.matrix.get(k, i) * v1);
						} else {
							w.set(0, i, w.get(0, i) + this.matrix.get(k, j) * this.matrix.get(k, i));
						}
					}
					
					//sets the final matrix
					this.matrix.set(j, i, this.matrix.get(j, i) - vect.get(0, j) * w.get(0, i));
					for (int k = j + 1; k < this.matrix.getRowDimension(); k++) {
						if (i > j) {
							this.matrix.set(k, i, this.matrix.get(k, i) - vect.get(0, j) * this.matrix.get(k, j) * w.get(0, i));
						}
					}
				}
			}
		}
		
		//code to return Q from the QR factorization
		public Matrix getQ() {
			int m = Math.max(matrix.getColumnDimension(), matrix.getRowDimension());
			Matrix q = new Matrix(m, m);
			
			for (int i = 0; i < Math.min(q.getRowDimension(), q.getColumnDimension()); i++) {
				q.set(i, i, 1);
			}
			
			for (int k = matrix.getColumnDimension() - 1; k >= 0; k--) {
				for ( int j = k; j < q.getColumnDimension(); j++) {
					
					double w = q.get(k, j);
					for (int i = k + 1; i < q.getRowDimension(); i++) {
						w += matrix.get(i, k) * q.get(i, j);
					}
					
					q.set(k, j, q.get(k, j) - vect.get(0, k) * w);
					for (int i = k + 1; i < q.getRowDimension(); i++) {
						q.set(i, j, q.get(i, j) - vect.get(0, k) * matrix.get(i, k) * w);
					}
				}
			}
			return q.getMatrix(0, q.getRowDimension() - 1, 0, 2);
		}
		
		//returns R from the QR factorization
		public Matrix getR() {
			Matrix r = matrix.copy();
			
			
			for (int i = 0; i < r.getRowDimension(); i++) {
				for (int j = 0; j < i && j < r.getColumnDimension(); j++) {
					r.set(i, j, 0);
				}
			}
			return r.getMatrix(0, 2, 0, 2);
		}
	}
