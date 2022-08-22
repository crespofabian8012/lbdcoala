package beast.evolution.substitutionmodel;

import java.util.Arrays;

import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.Nucleotide;
import beast.core.Description;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;

@Description("Jukes Cantor  for unphased diploid  genotypes")
public class JukesCantorGenotypes extends SubstitutionModel.Base {
	
	 public JukesCantorGenotypes() {
	        // this is added to avoid a parsing error inherited from superclass because frequencies are not provided.
	        frequenciesInput.setRule(Validate.OPTIONAL);
	        try {
	            // this call will be made twice when constructed from XML
	            // but this ensures that the object is validly constructed for testing purposes.
	            initAndValidate();
	        } catch (Exception e) {
	            e.printStackTrace();
	            throw new RuntimeException("initAndValidate() call failed when constructing JukesCantor()");
	        }
	    }
	
	 double[] frequencies;
	 EigenDecomposition eigenDecomposition;
	 
	 
	  @Override
	    public void initAndValidate() {
		    double minus_16_over_15 = - 16.0/ 15.0;
		    //eigen values
	        double[] eval = new double[]{0.0,                 -1.0666666666666666, -1.0666666666666666, -1.0666666666666666,
	        		                     -1.0666666666666666, -1.0666666666666666, -1.0666666666666666, -1.0666666666666666,
	        		                     -1.0666666666666666, -1.0666666666666666, -1.0666666666666666, -1.0666666666666666,
	        		                     -1.0666666666666666, -1.0666666666666666, -1.0666666666666666, -1.0666666666666666};
	        //eigen vectors
	        double[] evec = new double[]{1.0,  2.0,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                                         1.0, -2.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                                         1.0,  2.0, -0.5,  0.0,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                                         1.0, -2.0,  0.0, -0.5,  0.0,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                                         1.0,  2.0,  0.0,  0.0, -0.5,  0.0,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                                         1.0, -2.0,  0.0,  0.0,  0.0, -0.5,  0.0,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                                         1.0,  2.0,  0.0,  0.0,  0.0,  0.0, -0.5,  0.0,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                                         1.0, -2.0,  0.0,  0.0,  0.0,  0.0,  0.0, -0.5,  0.0,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                                         1.0,  2.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -0.5,  0.0,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0,
                                         1.0, -2.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -0.5,  0.0,  0.5,  0.0,  0.0,  0.0,  0.0,
                                         1.0,  2.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -0.5,  0.0,  0.5,  0.0,  0.0,  0.0,
                                         1.0, -2.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -0.5,  0.0,  0.5,  0.0,  0.0,
                                         1.0,  2.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -0.5,  0.0,  0.5,  0.0,
                                         1.0, -2.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -0.5,  0.0,  0.5,
                                         1.0,  2.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -0.5,  0.0,
                                         1.0, -2.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -0.5};
	        //inverse eigen vectors: the inverse of eigen vector matrix
	        double[] ivec = new double[]{ 0.0625,    0.0625,   0.0625,   0.0625,  0.0625,   0.0625,   0.0625,   0.0625,  0.0625,   0.0625,  0.0625,   0.0625,   0.0625,   0.0625,   0.0625,   0.0625,
	                                     0.03125,  -0.03125,  0.03125, -0.03125, 0.03125, -0.03125,  0.03125, -0.03125, 0.03125, -0.03125, 0.03125, -0.03125,  0.03125, -0.03125,  0.03125, -0.03125,
	                                        1.75,       0.0,    -0.25,      0.0,   -0.25,      0.0,    -0.25,      0.0,   -0.25,      0.0,   -0.25,      0.0,    -0.25,      0.0,    -0.25,      0.0,
	                                         0.0,      1.75,      0.0,    -0.25,     0.0,    -0.25,      0.0,    -0.25,     0.0,    -0.25,     0.0,    -0.25,      0.0,    -0.25,      0.0,    -0.25,
	                                         1.5,       0.0,      1.5,      0.0,    -0.5,      0.0,     -0.5,      0.0,    -0.5,      0.0,    -0.5,      0.0,     -0.5,      0.0,     -0.5,      0.0,
	                                         0.0,       1.5,      0.0,      1.5,     0.0,     -0.5,      0.0,     -0.5,     0.0,     -0.5,     0.0,     -0.5,      0.0,     -0.5,      0.0,     -0.5,
	                                        1.25,       0.0,     1.25,      0.0,    1.25,      0.0,    -0.75,      0.0,   -0.75,      0.0,   -0.75,      0.0,    -0.75,      0.0,    -0.75,      0.0,
	                                         0.0,      1.25,      0.0,     1.25,     0.0,     1.25,      0.0,    -0.75,     0.0,    -0.75,     0.0,    -0.75,      0.0,    -0.75,      0.0,    -0.75,
	                                         1.0,       0.0,      1.0,      0.0,     1.0,      0.0,      1.0,      0.0,    -1.0,      0.0,    -1.0,      0.0,     -1.0,      0.0,     -1.0,      0.0,
	                                         0.0,       1.0,      0.0,      1.0,     0.0,      1.0,      0.0,      1.0,     0.0,     -1.0,     0.0,     -1.0,      0.0,     -1.0,      0.0,     -1.0,
	                                        0.75,       0.0,     0.75,      0.0,    0.75,      0.0,     0.75,      0.0,    0.75,      0.0,   -1.25,      0.0,    -1.25,      0.0,    -1.25,      0.0,
	                                         0.0,      0.75,      0.0,     0.75,     0.0,     0.75,      0.0,     0.75,     0.0,     0.75,     0.0,    -1.25,      0.0,    -1.25,      0.0,    -1.25,
	                                         0.5,       0.0,      0.5,      0.0,     0.5,      0.0,      0.5,      0.0,     0.5,      0.0,     0.5,      0.0,     -1.5,      0.0,     -1.5,      0.0,
	                                         0.0,       0.5,      0.0,      0.5,     0.0,      0.5,      0.0,      0.5,     0.0,      0.5,     0.0,      0.5,      0.0,     -1.5,      0.0,     -1.5,
	                                        0.25,       0.0,     0.25,      0.0,    0.25,      0.0,     0.25,      0.0,    0.25,      0.0,    0.25,      0.0,     0.25,      0.0,    -1.75,      0.0,
	                                         0.0,       0.25,     0.0,     0.25,     0.0,     0.25,      0.0,     0.25,     0.0,     0.25,     0.0,     0.25,      0.0,     0.25,      0.0,    -1.75};

	        eigenDecomposition = new EigenDecomposition(evec, ivec, eval);

	        if (frequenciesInput.get() != null) {
	            throw new RuntimeException("Frequencies must not be specified in Jukes-Cantor16 model. They are assumed equal.");
	        }

	        frequencies = new double[]{0.0625, 0.0625, 0.0625, 0.0625, 
	        		                   0.0625, 0.0625, 0.0625, 0.0625, 
	        		                   0.0625, 0.0625, 0.0625, 0.0625,
	        		                   0.0625, 0.0625, 0.0625, 0.0625};
	    }

	  @Override
	    public double[] getFrequencies() {
	        return frequencies;
	    }

	@Override
	public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
		    double delta = 16.0 / 15.0 * (startTime - endTime);
	        double pStay = (1.0 + 15.0 * Math.exp(-delta * rate)) / 16.0;
	        double pMove = (1.0 - Math.exp(-delta * rate)) / 16.0;
	        // fill the matrix with move probabilities
	        Arrays.fill(matrix, pMove);
	        // fill the diagonal
	        for (int i = 0; i < 16; i++) {
	            matrix[i * 17] = pStay;
	        }

	}

    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {
        return eigenDecomposition;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof Nucleotide;
    }
}
