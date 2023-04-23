package lbdcoal;


import java.util.List;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.RealParameter;


/**
 * @author Fausto Fabian Crespo Fernandez
 */
@Description("Coalescent approximation to Limit Birth Death coalescent process")
@Citation(value="Crespo, F., Posada D., Carsten W. (2021).\n" +
        "Coalescent models derived from Birth Death processes.\n" +
        "Theoretical population biology).",
        year = 2021, firstAuthorSurname = "Crespo")
public class LimitBirthDeathApproximation extends PopulationFunction.Abstract {
	final public Input<RealParameter> scaledGrowthRateParameterInput = new Input<>("scaledGrowthRate",
			"growth rate scaled by present-day population size N0 ->+infinity. ");
	final public Input<RealParameter> timeOriginParameterInput = new Input<>("timeOrigin",
			"Time of origin in model time. ");
	final public Input<Double> KParameterInput = new Input<>("K",
			"parameter K in M_k approximation to Limit Birth-Death coalescent. ", 0.8);
	
	final public Input<RealParameter> scaledMutRateParameterInput = new Input<>("popSize",
			"mutation rate scaled by present-day population size N0 ->+infinity ");
	//
	// Public stuff
	//

	@Override
	public void initAndValidate() {
//	    	if (scaledGrowthRateParameterInput.get() <= 0.0) {
//	    		 throw new IllegalArgumentException("the scaled growth rate must be positive");
//	    	}

		if (timeOriginParameterInput.get() != null) {
			// TODO if not time of origin provided, draw from the conditional distribution
			timeOriginParameterInput.get().setBounds(0.0, 1000.0);
		}
//	        if (growthRateParameter.get() != null) {
//	            growthRateParameter.get().setBounds(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
//	        }
	}
	 /**
     * @return initial population size.
     */
    public double getTheta() {
        return scaledMutRateParameterInput.get().getValue();
    }



	/**
	 * @return scaled growth rate.
	 */
	public double getScaledGrowthRate() {
		gamma = scaledGrowthRateParameterInput.get().getValue();
		return gamma;
	}

	/**
	 * sets scaled growth rate .
	 *
	 * @param newScaledGrowthRate new scaled growth rate
	 */

	public void setScaledGrowthRate(double newScaledGrowthRate) {
		gamma = newScaledGrowthRate;
	}

	/**
	 * @return time of origin.
	 */
	public final double getT() {
		T = timeOriginParameterInput.get().getValue();
		return T;
	}

	/**
	 * sets time origin to r.
	 *
	 * @param r
	 */

	public void setT(double r) {
		T = r;
	}
    
	// Implementation of abstract methods
	@Override
	public double getPopSize(double t) {
		  
	        if (t == 0) {
	            return getTheta();
	        } else {
	            return getTheta() * (2.0 / getIntensity(t) );
	        }
	       // if( t > getT()) 
	       // {
	        	
	       // 	return Double.NEGATIVE_INFINITY;
	        	
	       // }
	        
	        
		
	}

	/**
	 * Calculates the integral 1/N(x) dx between start and finish.
	 */
	@Override
	public double getIntegral(double start, double finish) {
		double gamma = getScaledGrowthRate();
		double T = getT();
		//double K = getK();
		double a, b, expMinusGammaT, expGammaFinish;
		if (gamma >= 0.0) {
			expMinusGammaT = Math.exp(-gamma * T);
			expGammaFinish = Math.exp(gamma * finish);
			a = (K / gamma) * (1 - expMinusGammaT) - expMinusGammaT;
			b = 1 - (K / gamma) * (1 - expMinusGammaT);

			return (2.0 / K) * Math.log((a * expGammaFinish + b) / (1 - expGammaFinish * expMinusGammaT));
		} else {
			return Double.POSITIVE_INFINITY;
		}
	}

	@Override
	public double getIntensity(double t) {
		double gamma = getScaledGrowthRate();
		double T = getT();
		//double K = getK();
		double a, b;
		double expMinusGammaT, expGammat, oneMinusExpGammaT, intensity;
		if (gamma >= 0.0) {
			expMinusGammaT = Math.exp(-gamma * T);
			expGammat = Math.exp(gamma * t);
			oneMinusExpGammaT = 1 - expMinusGammaT;
			intensity = 2 * expGammat * Math.pow(oneMinusExpGammaT, 2) / Math.pow(1 - expMinusGammaT * expGammat, 2);
			intensity = intensity
					/ (1 + (K / gamma) * oneMinusExpGammaT * (expGammat - 1) / (1 - expMinusGammaT * expGammat));
			return intensity;
		} else {
			return Double.POSITIVE_INFINITY;
		}
	}

	@Override
	public double getInverseIntensity(double x) {

		double gamma = getScaledGrowthRate();
		double T = getT();
	    //double K = getK();
		double a, b, expKt, expMinusGammaT, expGammat, oneMinusExpGammaT, invItensity;
		if (gamma >= 0.0) {
			expMinusGammaT = Math.exp(-gamma * T);
			expGammat = Math.exp(gamma * x);
			oneMinusExpGammaT = 1 - expMinusGammaT;
			a = (K / gamma) * oneMinusExpGammaT - expMinusGammaT;
			b = 1 - (K / gamma) * oneMinusExpGammaT;
			expKt = Math.exp(K * x / 2.0);
			invItensity = 2 * expGammat * Math.pow(oneMinusExpGammaT, 2) / Math.pow(1 - expMinusGammaT * expGammat, 2);
			invItensity = (1.0 / gamma) * Math.log((expKt - b) / (expKt * expMinusGammaT + a));
			return invItensity;
		} else {
			return Double.POSITIVE_INFINITY;
		}
	}

	// Implementation of abstract methods

	@Override
	public List<String> getParameterIds() {
		List<String> paramIDs = new ArrayList<>();
		paramIDs.add(scaledGrowthRateParameterInput.get().getID());
		paramIDs.add(timeOriginParameterInput.get().getID());
		paramIDs.add(scaledMutRateParameterInput.get().getID());
		return paramIDs;

	}
	private double gamma;
	private double T;
	private double K;

	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
