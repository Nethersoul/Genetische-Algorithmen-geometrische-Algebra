/*
Gaigen 2.5 Test Suite
*/
/*
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/
using System;
namespace c3ga_ns {
public class TestSuite  :  c3ga
{ 
// Missing dependencies definitions:


/// <summary>Generates a random versor.
/// The scale is uniformly distributed over the interval [0 1).
/// The maximum non-zero grade-part is 'grade'.
/// 
/// Only the basis vectors marked in 'basisVectorBitmap' will be used to generate the versor/blade.
/// Use 'basisVectorBitmap = -1' (the default) to use all basisvectors.
/// </summary>
/// <returns>random_versor_dont_mangle_1_returns_mv_ex(arg1, scale, grade, basisVectorBitmap, 0.01, scale * 4.0)
/// </returns>
public static mv random_versor_dont_mangle_1_returns_mv(double scale, int grade, int basisVectorBitmap) {
	double minimumNorm = 0.01;
	double largestCoordinate = 4.0;
	return random_versor_dont_mangle_1_returns_mv_ex(scale, grade, basisVectorBitmap, minimumNorm, scale * largestCoordinate);
}

/// <summary>This version of random_versor_dont_mangle_1_returns_mv() has extra arguments which help to avoid
/// near-singular blades / versors.
/// 
/// Near-singular blades / versors are avoid by testing the norm and largest coordinate
/// of the random blade / versor. If the test does not pass, the function recursively
/// tries to generate another random blade / versor.
/// 
/// 'minimumNorm' is the minimum allowed norm of the blade/versor before scaling. 
/// 'minimumNorm' must be > 0.0 for versors.
/// 
/// 'largestCoordinate' is the largest coordinate allowed after scaling.
/// 
/// </summary>
/// <returns>random_versor_dont_mangle_1_returns_mv_ex(arg1, scale, grade, basisVectorBitmap, 0.01, scale * 4.0)
/// </returns>
public static mv random_versor_dont_mangle_1_returns_mv_ex(double scale, int _grade, int basisVectorBitmap, 
		double minimumNorm, double largestCoordinate) 
{
	mv randomVector = new mv();
	//, tmp1, tmp2;
	double[] randomValues = new double[3];
	//double n2, mul;
	int grade = _grade;
	
	// set IR (intermediate result) to 1
	mv IR = new mv (1.0);

	while (grade > 0) {	// loop until grade == 0
		// fill array with random values
		for (int i = 0; i < 3; i++) 
			randomValues[i] = ((basisVectorBitmap & (1 << i)) == 0)
				? 0.0 
				: (-1.0 + 2.0 * genrand());
		
		// make it a multivector:
		randomVector.Set(GroupBitmap.GRADE_1, randomValues);
		
		// multiply 
		IR = gp(IR, randomVector);
		
		// lower grade
		grade--;
	}
	
	// compute norm/multiplier, apply it, or recurse if we happened to create a near-null versor
	double n2 = norm_returns_scalar(IR);
	if ((double)Math.Abs(n2) > minimumNorm * minimumNorm) {
		if (n2 != 0.0) {
			double mul = scale * genrand() / n2;
			if (IR.LargestCoordinate() * mul < largestCoordinate)
				return gp(IR, mul);
		}
		else if (IR.LargestCoordinate() < largestCoordinate)
			return IR;
	}
	
	// try again:
	return random_versor_dont_mangle_1_returns_mv_ex(scale, _grade, basisVectorBitmap, minimumNorm, largestCoordinate); 
}
/// <summary>Returns random vector with a scale in the interval [0, scale)
/// </summary>
public static vector random_vector_dont_mangle_2_ex(double scale, double minimumNorm, double largestCoordinate)
{
	vector tmp = new vector();
	double n, mul, lc;
	bool nullBlade;
	double re1 = -1.0 + 2.0 * genrand(), re2 = -1.0 + 2.0 * genrand(), re3 = -1.0 + 2.0 * genrand();
	tmp.m_e1 = re1;
	tmp.m_e2 = re2;
	tmp.m_e3 = re3;

	n = norm_returns_scalar(tmp);
	lc = tmp.LargestCoordinate();
	nullBlade = ((n == 0.0) && (lc != 0.0));
	if ((n < minimumNorm) && (!nullBlade)) {
		return random_vector_dont_mangle_2_ex(scale, minimumNorm, largestCoordinate);
	}
	if (n < 0.0001) {
		mul = 1.0;
	}
	else {
		mul = scale * (-1.0 + 2.0 * genrand()) / n;
		if ((lc * Math.Abs(mul)) > largestCoordinate ) {
			return random_vector_dont_mangle_2_ex(scale, minimumNorm, largestCoordinate);
		}
	}
	return new vector(vector.coord_e1_e2_e3,
			mul*tmp.m_e1, // e1
			mul*tmp.m_e2, // e2
			mul*tmp.m_e3 // e3
		);
}
/// <summary>Returns random vector with a scale in the interval [0, scale)
/// </summary>
public static vector random_vector_dont_mangle_2(double scale)
{
	double minimumNorm = 0.0001;
	double largestCoordinate = 4.0;
	return random_vector_dont_mangle_2_ex(scale, minimumNorm, scale * largestCoordinate);
}
/// <summary>Returns random bivector with a scale in the interval [0, scale)
/// </summary>
public static bivector random_bivector_dont_mangle_4_ex(double scale, double minimumNorm, double largestCoordinate)
{
	bivector tmp = new bivector();
	double n, mul, lc;
	bool nullBlade;
	double re1_e2 = -1.0 + 2.0 * genrand(), re2_e3 = -1.0 + 2.0 * genrand(), re3_e1 = -1.0 + 2.0 * genrand();
	tmp.m_e1_e2 = re1_e2;
	tmp.m_e2_e3 = re2_e3;
	tmp.m_e3_e1 = re3_e1;

	n = norm_returns_scalar(tmp);
	lc = tmp.LargestCoordinate();
	nullBlade = ((n == 0.0) && (lc != 0.0));
	if ((n < minimumNorm) && (!nullBlade)) {
		return random_bivector_dont_mangle_4_ex(scale, minimumNorm, largestCoordinate);
	}
	if (n < 0.0001) {
		mul = 1.0;
	}
	else {
		mul = scale * (-1.0 + 2.0 * genrand()) / n;
		if ((lc * Math.Abs(mul)) > largestCoordinate ) {
			return random_bivector_dont_mangle_4_ex(scale, minimumNorm, largestCoordinate);
		}
	}
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			mul*tmp.m_e1_e2, // e1_e2
			mul*tmp.m_e2_e3, // e2_e3
			mul*tmp.m_e3_e1 // e3_e1
		);
}
/// <summary>Returns random bivector with a scale in the interval [0, scale)
/// </summary>
public static bivector random_bivector_dont_mangle_4(double scale)
{
	double minimumNorm = 0.0001;
	double largestCoordinate = 4.0;
	return random_bivector_dont_mangle_4_ex(scale, minimumNorm, scale * largestCoordinate);
}
/// <summary>Returns random trivector with a scale in the interval [0, scale)
/// </summary>
public static trivector random_trivector_dont_mangle_7_ex(double scale, double minimumNorm, double largestCoordinate)
{
	trivector tmp = new trivector();
	double n, mul, lc;
	bool nullBlade;
	double re1_e2_e3 = -1.0 + 2.0 * genrand();
	tmp.m_e1_e2_e3 = re1_e2_e3;

	n = norm_returns_scalar(tmp);
	lc = tmp.LargestCoordinate();
	nullBlade = ((n == 0.0) && (lc != 0.0));
	if ((n < minimumNorm) && (!nullBlade)) {
		return random_trivector_dont_mangle_7_ex(scale, minimumNorm, largestCoordinate);
	}
	if (n < 0.0001) {
		mul = 1.0;
	}
	else {
		mul = scale * (-1.0 + 2.0 * genrand()) / n;
		if ((lc * Math.Abs(mul)) > largestCoordinate ) {
			return random_trivector_dont_mangle_7_ex(scale, minimumNorm, largestCoordinate);
		}
	}
	return new trivector(trivector.coord_e1e2e3,
			mul*tmp.m_e1_e2_e3 // e1_e2_e3
		);
}
/// <summary>Returns random trivector with a scale in the interval [0, scale)
/// </summary>
public static trivector random_trivector_dont_mangle_7(double scale)
{
	double minimumNorm = 0.0001;
	double largestCoordinate = 4.0;
	return random_trivector_dont_mangle_7_ex(scale, minimumNorm, scale * largestCoordinate);
}
/// <summary>Returns random rotor with a scale in the interval [0, scale)
/// </summary>
public static rotor random_rotor_dont_mangle_8_ex(double scale, double minimumNorm, double largestCoordinate)
{
	rotor tmp = new rotor();
	double n, mul, lc;
	bool nullBlade;
	double r0_e1 = -1.0 + 2.0 * genrand(), r0_e2 = -1.0 + 2.0 * genrand(), r0_e3 = -1.0 + 2.0 * genrand(), 
			r1_e1 = -1.0 + 2.0 * genrand(), r1_e2 = -1.0 + 2.0 * genrand(), r1_e3 = -1.0 + 2.0 * genrand();
	tmp.m_scalar = (r0_e1*r1_e1+r0_e2*r1_e2+r0_e3*r1_e3);
	tmp.m_e1_e2 = (r0_e1*r1_e2-r0_e2*r1_e1);
	tmp.m_e2_e3 = (r0_e2*r1_e3-r0_e3*r1_e2);
	tmp.m_e3_e1 = -(r0_e1*r1_e3-r0_e3*r1_e1);

	n = norm_returns_scalar(tmp);
	lc = tmp.LargestCoordinate();
	nullBlade = false;
	if ((n < minimumNorm) && (!nullBlade)) {
		return random_rotor_dont_mangle_8_ex(scale, minimumNorm, largestCoordinate);
	}
	if (n < 0.0001) {
		mul = 1.0;
	}
	else {
		mul = scale * (-1.0 + 2.0 * genrand()) / n;
		if ((lc * Math.Abs(mul)) > largestCoordinate ) {
			return random_rotor_dont_mangle_8_ex(scale, minimumNorm, largestCoordinate);
		}
	}
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			mul*tmp.m_scalar, // scalar
			mul*tmp.m_e1_e2, // e1_e2
			mul*tmp.m_e2_e3, // e2_e3
			mul*tmp.m_e3_e1 // e3_e1
		);
}
/// <summary>Returns random rotor with a scale in the interval [0, scale)
/// </summary>
public static rotor random_rotor_dont_mangle_8(double scale)
{
	double minimumNorm = 0.0001;
	double largestCoordinate = 4.0;
	return random_rotor_dont_mangle_8_ex(scale, minimumNorm, scale * largestCoordinate);
}
/// <summary>Returns random e1_t with a scale in the interval [0, scale)
/// </summary>
public static e1_t random_e1_t_dont_mangle_10_ex(double scale, double minimumNorm, double largestCoordinate)
{
	return new e1_t(		);
}
/// <summary>Returns random e1_t with a scale in the interval [0, scale)
/// </summary>
public static e1_t random_e1_t_dont_mangle_10(double scale)
{
	double minimumNorm = 0.0001;
	double largestCoordinate = 4.0;
	return random_e1_t_dont_mangle_10_ex(scale, minimumNorm, scale * largestCoordinate);
}
/// <summary>Returns random e2_t with a scale in the interval [0, scale)
/// </summary>
public static e2_t random_e2_t_dont_mangle_11_ex(double scale, double minimumNorm, double largestCoordinate)
{
	return new e2_t(		);
}
/// <summary>Returns random e2_t with a scale in the interval [0, scale)
/// </summary>
public static e2_t random_e2_t_dont_mangle_11(double scale)
{
	double minimumNorm = 0.0001;
	double largestCoordinate = 4.0;
	return random_e2_t_dont_mangle_11_ex(scale, minimumNorm, scale * largestCoordinate);
}
/// <summary>Returns random I3_t with a scale in the interval [0, scale)
/// </summary>
public static I3_t random_I3_t_dont_mangle_13_ex(double scale, double minimumNorm, double largestCoordinate)
{
	return new I3_t(		);
}
/// <summary>Returns random I3_t with a scale in the interval [0, scale)
/// </summary>
public static I3_t random_I3_t_dont_mangle_13(double scale)
{
	double minimumNorm = 0.0001;
	double largestCoordinate = 4.0;
	return random_I3_t_dont_mangle_13_ex(scale, minimumNorm, scale * largestCoordinate);
}


/// <summary>Generates a random blade.
/// The scale is uniformly distributed over the interval [0 1).
/// The maximum non-zero grade-part is 'grade'.
/// 
/// Only the basis vectors marked in 'basisVectorBitmap' will be used to generate the versor/blade.
/// Use 'basisVectorBitmap = -1' (the default) to use all basisvectors.
/// </summary>
/// <returns>random_blade_dont_mangle_23_returns_mv_ex(arg1, scale, grade, basisVectorBitmap, 0.01, scale * 4.0)
/// </returns>
public static mv random_blade_dont_mangle_23_returns_mv(double scale, int grade, int basisVectorBitmap) {
	double minimumNorm = 0.01;
	double largestCoordinate = 4.0;
	return random_blade_dont_mangle_23_returns_mv_ex(scale, grade, basisVectorBitmap, minimumNorm, scale * largestCoordinate);
}

/// <summary>This version of random_blade_dont_mangle_23_returns_mv() has extra arguments which help to avoid
/// near-singular blades / versors.
/// 
/// Near-singular blades / versors are avoid by testing the norm and largest coordinate
/// of the random blade / versor. If the test does not pass, the function recursively
/// tries to generate another random blade / versor.
/// 
/// 'minimumNorm' is the minimum allowed norm of the blade/versor before scaling. 
/// 'minimumNorm' must be > 0.0 for versors.
/// 
/// 'largestCoordinate' is the largest coordinate allowed after scaling.
/// 
/// </summary>
/// <returns>random_blade_dont_mangle_23_returns_mv_ex(arg1, scale, grade, basisVectorBitmap, 0.01, scale * 4.0)
/// </returns>
public static mv random_blade_dont_mangle_23_returns_mv_ex(double scale, int _grade, int basisVectorBitmap, 
		double minimumNorm, double largestCoordinate) 
{
	mv randomVector = new mv();
	//, tmp1, tmp2;
	double[] randomValues = new double[3];
	//double n2, mul;
	int grade = _grade;
	
	// set IR (intermediate result) to 1
	mv IR = new mv (1.0);

	while (grade > 0) {	// loop until grade == 0
		// fill array with random values
		for (int i = 0; i < 3; i++) 
			randomValues[i] = ((basisVectorBitmap & (1 << i)) == 0)
				? 0.0 
				: (-1.0 + 2.0 * genrand());
		
		// make it a multivector:
		randomVector.Set(GroupBitmap.GRADE_1, randomValues);
		
		// multiply 
		IR = op(IR, randomVector);
		
		// lower grade
		grade--;
	}
	
	// compute norm/multiplier, apply it, or recurse if we happened to create a near-null versor
	double n2 = norm_returns_scalar(IR);
	if ((double)Math.Abs(n2) > minimumNorm * minimumNorm) {
		if (n2 != 0.0) {
			double mul = scale * genrand() / n2;
			if (IR.LargestCoordinate() * mul < largestCoordinate)
				return gp(IR, mul);
		}
		else if (IR.LargestCoordinate() < largestCoordinate)
			return IR;
	}
	
	// try again:
	return random_blade_dont_mangle_23_returns_mv_ex(scale, _grade, basisVectorBitmap, minimumNorm, largestCoordinate); 
}
/// <summary>Returns random e3_t with a scale in the interval [0, scale)
/// </summary>
public static e3_t random_e3_t_dont_mangle_75_ex(double scale, double minimumNorm, double largestCoordinate)
{
	return new e3_t(		);
}
/// <summary>Returns random e3_t with a scale in the interval [0, scale)
/// </summary>
public static e3_t random_e3_t_dont_mangle_75(double scale)
{
	double minimumNorm = 0.0001;
	double largestCoordinate = 4.0;
	return random_e3_t_dont_mangle_75_ex(scale, minimumNorm, scale * largestCoordinate);
}
/// <summary>Returns random oddVersor with a scale in the interval [0, scale)
/// </summary>
public static oddVersor random_oddVersor_dont_mangle_93_ex(double scale, double minimumNorm, double largestCoordinate)
{
	oddVersor tmp = new oddVersor();
	double n, mul, lc;
	bool nullBlade;
	double r0_e1 = -1.0 + 2.0 * genrand(), r0_e2 = -1.0 + 2.0 * genrand(), r0_e3 = -1.0 + 2.0 * genrand(), 
			r1_e1 = -1.0 + 2.0 * genrand(), r1_e2 = -1.0 + 2.0 * genrand(), r1_e3 = -1.0 + 2.0 * genrand(), 
			r2_e1 = -1.0 + 2.0 * genrand(), r2_e2 = -1.0 + 2.0 * genrand(), r2_e3 = -1.0 + 2.0 * genrand();
	tmp.m_e1 = (r0_e1*r1_e1*r2_e1+r0_e1*r1_e2*r2_e2+r0_e1*r1_e3*r2_e3-r0_e2*r1_e1*r2_e2+r0_e2*r1_e2*r2_e1-r0_e3*r1_e1*r2_e3+r0_e3*r1_e3*r2_e1);
	tmp.m_e2 = (r0_e1*r1_e1*r2_e2-r0_e1*r1_e2*r2_e1+r0_e2*r1_e1*r2_e1+r0_e2*r1_e2*r2_e2+r0_e2*r1_e3*r2_e3-r0_e3*r1_e2*r2_e3+r0_e3*r1_e3*r2_e2);
	tmp.m_e3 = (r0_e1*r1_e1*r2_e3-r0_e1*r1_e3*r2_e1+r0_e2*r1_e2*r2_e3-r0_e2*r1_e3*r2_e2+r0_e3*r1_e1*r2_e1+r0_e3*r1_e2*r2_e2+r0_e3*r1_e3*r2_e3);
	tmp.m_e1_e2_e3 = (r0_e1*r1_e2*r2_e3-r0_e1*r1_e3*r2_e2-r0_e2*r1_e1*r2_e3+r0_e2*r1_e3*r2_e1+r0_e3*r1_e1*r2_e2-r0_e3*r1_e2*r2_e1);

	n = norm_dont_mangle_381_returns_scalar(tmp);
	lc = tmp.LargestCoordinate();
	nullBlade = false;
	if ((n < minimumNorm) && (!nullBlade)) {
		return random_oddVersor_dont_mangle_93_ex(scale, minimumNorm, largestCoordinate);
	}
	if (n < 0.0001) {
		mul = 1.0;
	}
	else {
		mul = scale * (-1.0 + 2.0 * genrand()) / n;
		if ((lc * Math.Abs(mul)) > largestCoordinate ) {
			return random_oddVersor_dont_mangle_93_ex(scale, minimumNorm, largestCoordinate);
		}
	}
	return new oddVersor(oddVersor.coord_e1_e2_e3_e1e2e3,
			mul*tmp.m_e1, // e1
			mul*tmp.m_e2, // e2
			mul*tmp.m_e3, // e3
			mul*tmp.m_e1_e2_e3 // e1_e2_e3
		);
}
/// <summary>Returns random oddVersor with a scale in the interval [0, scale)
/// </summary>
public static oddVersor random_oddVersor_dont_mangle_93(double scale)
{
	double minimumNorm = 0.0001;
	double largestCoordinate = 4.0;
	return random_oddVersor_dont_mangle_93_ex(scale, minimumNorm, scale * largestCoordinate);
}
/// <summary>Returns norm of e2_t using default metric.
/// </summary>
public static double norm_dont_mangle_380(e2_t a)
{
	return Math.Abs(1.0);

}
/// <summary>internal conversion function
/// </summary>
public static double norm_dont_mangle_380_returns_scalar(e2_t a) {
	return norm_dont_mangle_380(a);
}
/// <summary>Returns norm of oddVersor using default metric.
/// </summary>
public static double norm_dont_mangle_381(oddVersor a)
{
	return Math.Abs(Math.Sqrt((a.m_e1*a.m_e1+a.m_e1_e2_e3*a.m_e1_e2_e3+a.m_e2*a.m_e2+a.m_e3*a.m_e3)));

}
/// <summary>internal conversion function
/// </summary>
public static double norm_dont_mangle_381_returns_scalar(oddVersor a) {
	return norm_dont_mangle_381(a);
}
// Testing code declarations:
// Testing code inline definitions:
// Testing code definitions:

static int test_metric_default_mv(int NB_TESTS_SCALER) 
{
	int i, j;
	double[] arr = new double[3];
	double dif;
	mv A;
	mv[] bv = new mv[3];
	double[] M = new double[]{
		1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0	}; // metric matrix

	// get all basis vectors

	c3ga.Zero_3(arr);
	arr[0] = 1.0;
	bv[0] = new mv(GroupBitmap.GROUP_1, arr);

	c3ga.Zero_3(arr);
	arr[1] = 1.0;
	bv[1] = new mv(GroupBitmap.GROUP_1, arr);

	c3ga.Zero_3(arr);
	arr[2] = 1.0;
	bv[2] = new mv(GroupBitmap.GROUP_1, arr);

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			A = gp(bv[i], bv[j]);
			dif = M[i * 3 + j] - A.get_scalar();
			if ((dif < -1E-14) || (dif > 1E-14)) {
				Console.WriteLine("test_metric_default_mv() test failed for " + BasisVectorNames[i] + " " + BasisVectorNames[j]);
				return 0;
			}
		}
	}
	
	return 1;
}

static int test_parse_mv(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 32;
	mv A, B, C;
	String str;
	
	int i, basisVectorBitmap = -1;

	for (i = 0; i < NB_LOOPS; i++) {
		A = random_versor_dont_mangle_1_returns_mv(genrand(), (int)(genrand() * 3.5), basisVectorBitmap);
		
		str = A.ToString(
			"E20"
		
		);
		
		try {
			B = Parse(str);
		} catch (ParseException ex) {
			Console.WriteLine("Parse() test failed: " + ex.ToString());
			return 0; // failure
		}

		C = subtract(A, B);

		if (C.LargestCoordinate() > 1E-14) {
			Console.WriteLine("Parse() test failed (largest coordinate: " + C.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	
	return 1; // success
}

static int test_genrand_double(int NB_TESTS_SCALER) 
{
	int NB_BINS = 256;
	int NB_LOOPS = NB_BINS * 1024;
	int[] bins = new int[256];
	double avg = 0.0;
	double r;
	int idx, i;
	
	// generate a lot of random values, compute average (should be 0.5) and fill bins
	for (i = 0; i < NB_LOOPS; i++) {
		r = genrand();
		avg += r;
		idx = (int)(r * (double)NB_BINS);
		if (idx >= NB_BINS) idx = NB_BINS-1;
		bins[idx]++;
	}
	avg /= (double)NB_LOOPS;
	
	if ((avg < 0.49) || (avg > 0.51)) {
		Console.WriteLine("Random number generator genrand() likely failed: average is " + (double)avg + " instead of 0.5.");
		return 0; // failure
	}
	
	for (i = 0; i < NB_BINS; i++) {
		if ((bins[i] < (6 * NB_LOOPS / NB_BINS / 8)) ||
			(bins[i] > (10 * NB_LOOPS / NB_BINS / 8))) {
			Console.WriteLine("It is very likely that the random number generator genrand() failed:");
			Console.WriteLine("The distribution is not uniform (bin " + i + " of " + NB_BINS + ", value=" + bins[i] + ", expected value=" + (NB_LOOPS / NB_BINS) + ")");
			return 0; // failure
		}
	}
	
	return 1; // success
}

static int test_add_dont_mangle_382(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, C;
	int i, g;
	double s;
	int basisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		g = (int)(genrand() * 3.5);
		A = random_versor_dont_mangle_1_returns_mv(s, g, basisVectorBitmap);
		
		// B = -A
		B = negate(A);
		
		C = add(A, B);
		
		// check
		if (C.LargestCoordinate() > 1E-13) {
			Console.WriteLine("add() test failed");
			return 0; // failure
		}
		
	}
	return 1; // success
}

static int test_add_dont_mangle_387(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 6;
	vector A;
	vector B;
	vector C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_vector_dont_mangle_2(s);
		B = random_vector_dont_mangle_2(s);
		
		// add/subtract A and B
		
		C = add(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = add(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("add() test failed");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_add_dont_mangle_385(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 6;
	bivector A;
	bivector B;
	bivector C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_bivector_dont_mangle_4(s);
		B = random_bivector_dont_mangle_4(s);
		
		// add/subtract A and B
		
		C = add(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = add(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("add() test failed");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_add_dont_mangle_388(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	vector A;
	trivector B;
	oddVersor C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_vector_dont_mangle_2(s);
		B = random_trivector_dont_mangle_7(s);
		
		// add/subtract A and B
		
		C = add(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = add(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("add() test failed");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_add_dont_mangle_384(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 7;
	rotor A;
	bivector B;
	rotor C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_rotor_dont_mangle_8(s);
		B = random_bivector_dont_mangle_4(s);
		
		// add/subtract A and B
		
		C = add(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = add(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("add() test failed");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_add_dont_mangle_386(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 2;
	e1_t A;
	e2_t B;
	vector C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_e1_t_dont_mangle_10(s);
		B = random_e2_t_dont_mangle_11(s);
		
		// add/subtract A and B
		
		C = add(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = add(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("add() test failed");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_add_dont_mangle_383(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 2;
	e1_t A;
	I3_t B;
	oddVersor C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_e1_t_dont_mangle_10(s);
		B = random_I3_t_dont_mangle_13(s);
		
		// add/subtract A and B
		
		C = add(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = add(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("add() test failed");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_subtract_dont_mangle_389(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, C;
	int i, g;
	double s;
	int basisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		g = (int)(genrand() * 3.5);
		A = random_versor_dont_mangle_1_returns_mv(s, g, basisVectorBitmap);
		
		B = A;
		
		C = subtract(A, B);
		
		// check
		if (C.LargestCoordinate() > 1E-13) {
			Console.WriteLine("subtract() test failed");
			return 0; // failure
		}
		
	}
	return 1; // success
}

static int test_subtract_dont_mangle_390(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 6;
	vector A;
	vector B;
	vector C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_vector_dont_mangle_2(s);
		B = random_vector_dont_mangle_2(s);
		
		// add/subtract A and B
		
		C = subtract(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = subtract(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("subtract() test failed");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_subtract_dont_mangle_391(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 6;
	bivector A;
	bivector B;
	bivector C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_bivector_dont_mangle_4(s);
		B = random_bivector_dont_mangle_4(s);
		
		// add/subtract A and B
		
		C = subtract(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = subtract(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("subtract() test failed");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_subtract_dont_mangle_392(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 7;
	bivector A;
	rotor B;
	rotor C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_bivector_dont_mangle_4(s);
		B = random_rotor_dont_mangle_8(s);
		
		// add/subtract A and B
		
		C = subtract(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = subtract(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("subtract() test failed");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_subtract_dont_mangle_393(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	vector A;
	trivector B;
	oddVersor C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_vector_dont_mangle_2(s);
		B = random_trivector_dont_mangle_7(s);
		
		// add/subtract A and B
		
		C = subtract(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = subtract(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("subtract() test failed");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_applyVersor_dont_mangle_394(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	double baseScale = 0.5;
	int g, i;
	double s;
	mv V, IV, X, Y, VX, VY, IVVX, XdY, VXdVY, dif;
	int versorBasisVectorBitmap = 7; // note: random versors restricted to Euclidean basis vectors.
	int bladeBasisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor and its inverse. Optionally make versor unit.
		s = baseScale + genrand();
		g = (int)(genrand() * 3.5);
		V = random_versor_dont_mangle_1_returns_mv(s, g, versorBasisVectorBitmap);
		
		// avoid 'near'-singular versors
		if (V.LargestCoordinate() > 2.0)
			continue;		
		
		IV = versorInverse(V);

		// get two random blades		
		s = baseScale + genrand();
		g = (int)(genrand() * 3.5);
		X = random_blade_dont_mangle_23_returns_mv(s, g, bladeBasisVectorBitmap);
		s = baseScale + genrand();
		Y = random_blade_dont_mangle_23_returns_mv(s, g, bladeBasisVectorBitmap);

		// apply versor to blades
		VX = new mv(applyVersor(V, X));
		VY = new mv(applyVersor(V, Y));
		
		// compute inner product
		XdY = mhip(X, Y);
		VXdVY = mhip(VX, VY);
		
		// see if inner products are equal (versor application should not change metric)
		dif = subtract(XdY, VXdVY);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyVersor() test failed (metric test) (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
		
		// apply inverse transformation to VX
		IVVX = applyVersor(IV, VX);
		
		// see if X equals IVVX
		dif = subtract(X, IVVX);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyVersor() test failed (inverse test) (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_applyUnitVersor_dont_mangle_395(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	double baseScale = 0.5;
	int g, i;
	double s;
	mv V, IV, X, Y, VX, VY, IVVX, XdY, VXdVY, dif;
	mv tmp;
	int versorBasisVectorBitmap = 7; // note: random versors restricted to Euclidean basis vectors.
	int bladeBasisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor and its inverse. Optionally make versor unit.
		s = baseScale + genrand();
		g = (int)(genrand() * 3.5);
		V = random_versor_dont_mangle_1_returns_mv(s, g, versorBasisVectorBitmap);
		
		// avoid 'near'-singular versors
		if (V.LargestCoordinate() > 2.0)
			continue;		
		
		tmp = unit(V);
		V = tmp;
		IV = versorInverse(V);

		// get two random blades		
		s = baseScale + genrand();
		g = (int)(genrand() * 3.5);
		X = random_blade_dont_mangle_23_returns_mv(s, g, bladeBasisVectorBitmap);
		s = baseScale + genrand();
		Y = random_blade_dont_mangle_23_returns_mv(s, g, bladeBasisVectorBitmap);

		// apply versor to blades
		VX = new mv(applyUnitVersor(V, X));
		VY = new mv(applyUnitVersor(V, Y));
		
		// compute inner product
		XdY = mhip(X, Y);
		VXdVY = mhip(VX, VY);
		
		// see if inner products are equal (versor application should not change metric)
		dif = subtract(XdY, VXdVY);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyUnitVersor() test failed (metric test) (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
		
		// apply inverse transformation to VX
		IVVX = applyUnitVersor(IV, VX);
		
		// see if X equals IVVX
		dif = subtract(X, IVVX);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyUnitVersor() test failed (inverse test) (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_applyVersorWI_dont_mangle_396(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	double baseScale = 0.5;
	int g, i;
	double s;
	mv V, IV, X, Y, VX, VY, IVVX, XdY, VXdVY, dif;
	mv tmp;
	int versorBasisVectorBitmap = 7; // note: random versors restricted to Euclidean basis vectors.
	int bladeBasisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor and its inverse. Optionally make versor unit.
		s = baseScale + genrand();
		g = (int)(genrand() * 3.5);
		V = random_versor_dont_mangle_1_returns_mv(s, g, versorBasisVectorBitmap);
		
		// avoid 'near'-singular versors
		if (V.LargestCoordinate() > 2.0)
			continue;		
		
		tmp = unit(V);
		V = tmp;
		IV = versorInverse(V);

		// get two random blades		
		s = baseScale + genrand();
		g = (int)(genrand() * 3.5);
		X = random_blade_dont_mangle_23_returns_mv(s, g, bladeBasisVectorBitmap);
		s = baseScale + genrand();
		Y = random_blade_dont_mangle_23_returns_mv(s, g, bladeBasisVectorBitmap);

		// apply versor to blades
		VX = new mv(applyVersorWI(V, X, IV));
		VY = new mv(applyVersorWI(V, Y, IV));
		
		// compute inner product
		XdY = mhip(X, Y);
		VXdVY = mhip(VX, VY);
		
		// see if inner products are equal (versor application should not change metric)
		dif = subtract(XdY, VXdVY);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyVersorWI() test failed (metric test) (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
		
		// apply inverse transformation to VX
		IVVX = applyVersorWI(IV, VX, V);
		
		// see if X equals IVVX
		dif = subtract(X, IVVX);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyVersorWI() test failed (inverse test) (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_applyVersor_dont_mangle_397(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor V;
	rotor IV;
	vector X;
	vector VX, tmpVX2;
	mv gmvV, gmvX, gmvVX, gmvVX2, dif; // gmvIV, 
	int i;

	for (i = 0; i < NB_LOOPS; i++) {
		// generate random versor of target type, make unit if required, invert
		V = random_rotor_dont_mangle_8(genrand());
		
		IV = versorInverse(V);

		// avoid near-singular versors
		if ((V.LargestCoordinate() > 2.0) ||
			(IV.LargestCoordinate() > 2.0))
			continue;
		 

		// generate random SMV 
		X = random_vector_dont_mangle_2(genrand());

		//  apply random versor to random SMV, convert to GMV
		VX = applyVersor(V, X);
		gmvVX2 = new mv(VX);

		// convert all to GMV type, apply versor too as GMV
		gmvV = new mv(V);
//		gmvIV = new mv(IV);
		gmvX = new mv(X);
		gmvVX = applyVersor(gmvV, gmvX);
		
		// convert GMV back and forth to return type to fix possible constant coordinates
		tmpVX2 = new vector(gmvVX);
		gmvVX = new mv(tmpVX2);
		
		// see if VX equals gmvVX
		dif = subtract(gmvVX, gmvVX2);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyVersor() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}


static int test_applyUnitVersor_dont_mangle_398(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor V;
	rotor tmp;
	rotor IV;
	vector X;
	vector VX, tmpVX2;
	mv gmvV, gmvX, gmvVX, gmvVX2, dif; // gmvIV, 
	int i;

	for (i = 0; i < NB_LOOPS; i++) {
		// generate random versor of target type, make unit if required, invert
		V = random_rotor_dont_mangle_8(genrand());
		
		// Test if norm2(V) is positive, otherwise do not perform the test.
		// (because with a negative norm2(v), the reverse is not the inverse)
		if (norm2_returns_scalar(V) <= 0.0) continue;

		tmp = unit(V);
		V = tmp;
		IV = versorInverse(V);

		// avoid near-singular versors
		if ((V.LargestCoordinate() > 2.0) ||
			(IV.LargestCoordinate() > 2.0))
			continue;
		 

		// generate random SMV 
		X = random_vector_dont_mangle_2(genrand());

		//  apply random versor to random SMV, convert to GMV
		VX = applyUnitVersor(V, X);
		gmvVX2 = new mv(VX);

		// convert all to GMV type, apply versor too as GMV
		gmvV = new mv(V);
//		gmvIV = new mv(IV);
		gmvX = new mv(X);
		gmvVX = applyVersor(gmvV, gmvX);
		
		// convert GMV back and forth to return type to fix possible constant coordinates
		tmpVX2 = new vector(gmvVX);
		gmvVX = new mv(tmpVX2);
		
		// see if VX equals gmvVX
		dif = subtract(gmvVX, gmvVX2);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyUnitVersor() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}


static int test_applyVersorWI_dont_mangle_399(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor V;
	rotor tmp;
	rotor IV;
	vector X;
	vector VX, tmpVX2;
	mv gmvV, gmvX, gmvVX, gmvVX2, dif; // gmvIV, 
	int i;

	for (i = 0; i < NB_LOOPS; i++) {
		// generate random versor of target type, make unit if required, invert
		V = random_rotor_dont_mangle_8(genrand());
		
		// Test if norm2(V) is positive, otherwise do not perform the test.
		// (because with a negative norm2(v), the reverse is not the inverse)
		if (norm2_returns_scalar(V) <= 0.0) continue;

		tmp = unit(V);
		V = tmp;
		IV = versorInverse(V);

		// avoid near-singular versors
		if ((V.LargestCoordinate() > 2.0) ||
			(IV.LargestCoordinate() > 2.0))
			continue;
		 

		// generate random SMV 
		X = random_vector_dont_mangle_2(genrand());

		//  apply random versor to random SMV, convert to GMV
		VX = applyVersorWI(V, X, IV);
		gmvVX2 = new mv(VX);

		// convert all to GMV type, apply versor too as GMV
		gmvV = new mv(V);
//		gmvIV = new mv(IV);
		gmvX = new mv(X);
		gmvVX = applyVersor(gmvV, gmvX);
		
		// convert GMV back and forth to return type to fix possible constant coordinates
		tmpVX2 = new vector(gmvVX);
		gmvVX = new mv(tmpVX2);
		
		// see if VX equals gmvVX
		dif = subtract(gmvVX, gmvVX2);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyVersorWI() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}


static int test_applyVersor_dont_mangle_400(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor V;
	rotor IV;
	bivector X;
	bivector VX, tmpVX2;
	mv gmvV, gmvX, gmvVX, gmvVX2, dif; // gmvIV, 
	int i;

	for (i = 0; i < NB_LOOPS; i++) {
		// generate random versor of target type, make unit if required, invert
		V = random_rotor_dont_mangle_8(genrand());
		
		IV = versorInverse(V);

		// avoid near-singular versors
		if ((V.LargestCoordinate() > 2.0) ||
			(IV.LargestCoordinate() > 2.0))
			continue;
		 

		// generate random SMV 
		X = random_bivector_dont_mangle_4(genrand());

		//  apply random versor to random SMV, convert to GMV
		VX = applyVersor(V, X);
		gmvVX2 = new mv(VX);

		// convert all to GMV type, apply versor too as GMV
		gmvV = new mv(V);
//		gmvIV = new mv(IV);
		gmvX = new mv(X);
		gmvVX = applyVersor(gmvV, gmvX);
		
		// convert GMV back and forth to return type to fix possible constant coordinates
		tmpVX2 = new bivector(gmvVX);
		gmvVX = new mv(tmpVX2);
		
		// see if VX equals gmvVX
		dif = subtract(gmvVX, gmvVX2);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyVersor() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}


static int test_applyUnitVersor_dont_mangle_401(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor V;
	rotor tmp;
	rotor IV;
	bivector X;
	bivector VX, tmpVX2;
	mv gmvV, gmvX, gmvVX, gmvVX2, dif; // gmvIV, 
	int i;

	for (i = 0; i < NB_LOOPS; i++) {
		// generate random versor of target type, make unit if required, invert
		V = random_rotor_dont_mangle_8(genrand());
		
		// Test if norm2(V) is positive, otherwise do not perform the test.
		// (because with a negative norm2(v), the reverse is not the inverse)
		if (norm2_returns_scalar(V) <= 0.0) continue;

		tmp = unit(V);
		V = tmp;
		IV = versorInverse(V);

		// avoid near-singular versors
		if ((V.LargestCoordinate() > 2.0) ||
			(IV.LargestCoordinate() > 2.0))
			continue;
		 

		// generate random SMV 
		X = random_bivector_dont_mangle_4(genrand());

		//  apply random versor to random SMV, convert to GMV
		VX = applyUnitVersor(V, X);
		gmvVX2 = new mv(VX);

		// convert all to GMV type, apply versor too as GMV
		gmvV = new mv(V);
//		gmvIV = new mv(IV);
		gmvX = new mv(X);
		gmvVX = applyVersor(gmvV, gmvX);
		
		// convert GMV back and forth to return type to fix possible constant coordinates
		tmpVX2 = new bivector(gmvVX);
		gmvVX = new mv(tmpVX2);
		
		// see if VX equals gmvVX
		dif = subtract(gmvVX, gmvVX2);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyUnitVersor() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}


static int test_applyVersorWI_dont_mangle_402(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor V;
	rotor tmp;
	rotor IV;
	bivector X;
	bivector VX, tmpVX2;
	mv gmvV, gmvX, gmvVX, gmvVX2, dif; // gmvIV, 
	int i;

	for (i = 0; i < NB_LOOPS; i++) {
		// generate random versor of target type, make unit if required, invert
		V = random_rotor_dont_mangle_8(genrand());
		
		// Test if norm2(V) is positive, otherwise do not perform the test.
		// (because with a negative norm2(v), the reverse is not the inverse)
		if (norm2_returns_scalar(V) <= 0.0) continue;

		tmp = unit(V);
		V = tmp;
		IV = versorInverse(V);

		// avoid near-singular versors
		if ((V.LargestCoordinate() > 2.0) ||
			(IV.LargestCoordinate() > 2.0))
			continue;
		 

		// generate random SMV 
		X = random_bivector_dont_mangle_4(genrand());

		//  apply random versor to random SMV, convert to GMV
		VX = applyVersorWI(V, X, IV);
		gmvVX2 = new mv(VX);

		// convert all to GMV type, apply versor too as GMV
		gmvV = new mv(V);
//		gmvIV = new mv(IV);
		gmvX = new mv(X);
		gmvVX = applyVersor(gmvV, gmvX);
		
		// convert GMV back and forth to return type to fix possible constant coordinates
		tmpVX2 = new bivector(gmvVX);
		gmvVX = new mv(tmpVX2);
		
		// see if VX equals gmvVX
		dif = subtract(gmvVX, gmvVX2);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyVersorWI() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}


static int test_applyVersor_dont_mangle_403(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor V;
	rotor IV;
	trivector X;
	trivector VX, tmpVX2;
	mv gmvV, gmvX, gmvVX, gmvVX2, dif; // gmvIV, 
	int i;

	for (i = 0; i < NB_LOOPS; i++) {
		// generate random versor of target type, make unit if required, invert
		V = random_rotor_dont_mangle_8(genrand());
		
		IV = versorInverse(V);

		// avoid near-singular versors
		if ((V.LargestCoordinate() > 2.0) ||
			(IV.LargestCoordinate() > 2.0))
			continue;
		 

		// generate random SMV 
		X = random_trivector_dont_mangle_7(genrand());

		//  apply random versor to random SMV, convert to GMV
		VX = applyVersor(V, X);
		gmvVX2 = new mv(VX);

		// convert all to GMV type, apply versor too as GMV
		gmvV = new mv(V);
//		gmvIV = new mv(IV);
		gmvX = new mv(X);
		gmvVX = applyVersor(gmvV, gmvX);
		
		// convert GMV back and forth to return type to fix possible constant coordinates
		tmpVX2 = new trivector(gmvVX);
		gmvVX = new mv(tmpVX2);
		
		// see if VX equals gmvVX
		dif = subtract(gmvVX, gmvVX2);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyVersor() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}


static int test_applyUnitVersor_dont_mangle_404(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor V;
	rotor tmp;
	rotor IV;
	trivector X;
	trivector VX, tmpVX2;
	mv gmvV, gmvX, gmvVX, gmvVX2, dif; // gmvIV, 
	int i;

	for (i = 0; i < NB_LOOPS; i++) {
		// generate random versor of target type, make unit if required, invert
		V = random_rotor_dont_mangle_8(genrand());
		
		// Test if norm2(V) is positive, otherwise do not perform the test.
		// (because with a negative norm2(v), the reverse is not the inverse)
		if (norm2_returns_scalar(V) <= 0.0) continue;

		tmp = unit(V);
		V = tmp;
		IV = versorInverse(V);

		// avoid near-singular versors
		if ((V.LargestCoordinate() > 2.0) ||
			(IV.LargestCoordinate() > 2.0))
			continue;
		 

		// generate random SMV 
		X = random_trivector_dont_mangle_7(genrand());

		//  apply random versor to random SMV, convert to GMV
		VX = applyUnitVersor(V, X);
		gmvVX2 = new mv(VX);

		// convert all to GMV type, apply versor too as GMV
		gmvV = new mv(V);
//		gmvIV = new mv(IV);
		gmvX = new mv(X);
		gmvVX = applyVersor(gmvV, gmvX);
		
		// convert GMV back and forth to return type to fix possible constant coordinates
		tmpVX2 = new trivector(gmvVX);
		gmvVX = new mv(tmpVX2);
		
		// see if VX equals gmvVX
		dif = subtract(gmvVX, gmvVX2);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyUnitVersor() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}


static int test_applyVersorWI_dont_mangle_405(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor V;
	rotor tmp;
	rotor IV;
	trivector X;
	trivector VX, tmpVX2;
	mv gmvV, gmvX, gmvVX, gmvVX2, dif; // gmvIV, 
	int i;

	for (i = 0; i < NB_LOOPS; i++) {
		// generate random versor of target type, make unit if required, invert
		V = random_rotor_dont_mangle_8(genrand());
		
		// Test if norm2(V) is positive, otherwise do not perform the test.
		// (because with a negative norm2(v), the reverse is not the inverse)
		if (norm2_returns_scalar(V) <= 0.0) continue;

		tmp = unit(V);
		V = tmp;
		IV = versorInverse(V);

		// avoid near-singular versors
		if ((V.LargestCoordinate() > 2.0) ||
			(IV.LargestCoordinate() > 2.0))
			continue;
		 

		// generate random SMV 
		X = random_trivector_dont_mangle_7(genrand());

		//  apply random versor to random SMV, convert to GMV
		VX = applyVersorWI(V, X, IV);
		gmvVX2 = new mv(VX);

		// convert all to GMV type, apply versor too as GMV
		gmvV = new mv(V);
//		gmvIV = new mv(IV);
		gmvX = new mv(X);
		gmvVX = applyVersor(gmvV, gmvX);
		
		// convert GMV back and forth to return type to fix possible constant coordinates
		tmpVX2 = new trivector(gmvVX);
		gmvVX = new mv(tmpVX2);
		
		// see if VX equals gmvVX
		dif = subtract(gmvVX, gmvVX2);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyVersorWI() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}


static int test_applyVersor_dont_mangle_406(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor V;
	rotor IV;
	e1_t X;
	vector VX, tmpVX2;
	mv gmvV, gmvX, gmvVX, gmvVX2, dif; // gmvIV, 
	int i;

	for (i = 0; i < NB_LOOPS; i++) {
		// generate random versor of target type, make unit if required, invert
		V = random_rotor_dont_mangle_8(genrand());
		
		IV = versorInverse(V);

		// avoid near-singular versors
		if ((V.LargestCoordinate() > 2.0) ||
			(IV.LargestCoordinate() > 2.0))
			continue;
		 

		// generate random SMV 
		X = random_e1_t_dont_mangle_10(genrand());

		//  apply random versor to random SMV, convert to GMV
		VX = applyVersor(V, X);
		gmvVX2 = new mv(VX);

		// convert all to GMV type, apply versor too as GMV
		gmvV = new mv(V);
//		gmvIV = new mv(IV);
		gmvX = new mv(X);
		gmvVX = applyVersor(gmvV, gmvX);
		
		// convert GMV back and forth to return type to fix possible constant coordinates
		tmpVX2 = new vector(gmvVX);
		gmvVX = new mv(tmpVX2);
		
		// see if VX equals gmvVX
		dif = subtract(gmvVX, gmvVX2);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyVersor() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}


static int test_applyUnitVersor_dont_mangle_407(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor V;
	rotor tmp;
	rotor IV;
	e2_t X;
	vector VX, tmpVX2;
	mv gmvV, gmvX, gmvVX, gmvVX2, dif; // gmvIV, 
	int i;

	for (i = 0; i < NB_LOOPS; i++) {
		// generate random versor of target type, make unit if required, invert
		V = random_rotor_dont_mangle_8(genrand());
		
		// Test if norm2(V) is positive, otherwise do not perform the test.
		// (because with a negative norm2(v), the reverse is not the inverse)
		if (norm2_returns_scalar(V) <= 0.0) continue;

		tmp = unit(V);
		V = tmp;
		IV = versorInverse(V);

		// avoid near-singular versors
		if ((V.LargestCoordinate() > 2.0) ||
			(IV.LargestCoordinate() > 2.0))
			continue;
		 

		// generate random SMV 
		X = random_e2_t_dont_mangle_11(genrand());

		//  apply random versor to random SMV, convert to GMV
		VX = applyUnitVersor(V, X);
		gmvVX2 = new mv(VX);

		// convert all to GMV type, apply versor too as GMV
		gmvV = new mv(V);
//		gmvIV = new mv(IV);
		gmvX = new mv(X);
		gmvVX = applyVersor(gmvV, gmvX);
		
		// convert GMV back and forth to return type to fix possible constant coordinates
		tmpVX2 = new vector(gmvVX);
		gmvVX = new mv(tmpVX2);
		
		// see if VX equals gmvVX
		dif = subtract(gmvVX, gmvVX2);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyUnitVersor() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}


static int test_applyVersor_dont_mangle_408(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor V;
	rotor IV;
	I3_t X;
	trivector VX, tmpVX2;
	mv gmvV, gmvX, gmvVX, gmvVX2, dif; // gmvIV, 
	int i;

	for (i = 0; i < NB_LOOPS; i++) {
		// generate random versor of target type, make unit if required, invert
		V = random_rotor_dont_mangle_8(genrand());
		
		IV = versorInverse(V);

		// avoid near-singular versors
		if ((V.LargestCoordinate() > 2.0) ||
			(IV.LargestCoordinate() > 2.0))
			continue;
		 

		// generate random SMV 
		X = random_I3_t_dont_mangle_13(genrand());

		//  apply random versor to random SMV, convert to GMV
		VX = applyVersor(V, X);
		gmvVX2 = new mv(VX);

		// convert all to GMV type, apply versor too as GMV
		gmvV = new mv(V);
//		gmvIV = new mv(IV);
		gmvX = new mv(X);
		gmvVX = applyVersor(gmvV, gmvX);
		
		// convert GMV back and forth to return type to fix possible constant coordinates
		tmpVX2 = new trivector(gmvVX);
		gmvVX = new mv(tmpVX2);
		
		// see if VX equals gmvVX
		dif = subtract(gmvVX, gmvVX2);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyVersor() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}


static int test_applyUnitVersor_dont_mangle_409(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor V;
	rotor tmp;
	rotor IV;
	I3_t X;
	trivector VX, tmpVX2;
	mv gmvV, gmvX, gmvVX, gmvVX2, dif; // gmvIV, 
	int i;

	for (i = 0; i < NB_LOOPS; i++) {
		// generate random versor of target type, make unit if required, invert
		V = random_rotor_dont_mangle_8(genrand());
		
		// Test if norm2(V) is positive, otherwise do not perform the test.
		// (because with a negative norm2(v), the reverse is not the inverse)
		if (norm2_returns_scalar(V) <= 0.0) continue;

		tmp = unit(V);
		V = tmp;
		IV = versorInverse(V);

		// avoid near-singular versors
		if ((V.LargestCoordinate() > 2.0) ||
			(IV.LargestCoordinate() > 2.0))
			continue;
		 

		// generate random SMV 
		X = random_I3_t_dont_mangle_13(genrand());

		//  apply random versor to random SMV, convert to GMV
		VX = applyUnitVersor(V, X);
		gmvVX2 = new mv(VX);

		// convert all to GMV type, apply versor too as GMV
		gmvV = new mv(V);
//		gmvIV = new mv(IV);
		gmvX = new mv(X);
		gmvVX = applyVersor(gmvV, gmvX);
		
		// convert GMV back and forth to return type to fix possible constant coordinates
		tmpVX2 = new trivector(gmvVX);
		gmvVX = new mv(tmpVX2);
		
		// see if VX equals gmvVX
		dif = subtract(gmvVX, gmvVX2);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyUnitVersor() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}


static int test_applyVersorWI_dont_mangle_410(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor V;
	rotor tmp;
	rotor IV;
	I3_t X;
	trivector VX, tmpVX2;
	mv gmvV, gmvX, gmvVX, gmvVX2, dif; // gmvIV, 
	int i;

	for (i = 0; i < NB_LOOPS; i++) {
		// generate random versor of target type, make unit if required, invert
		V = random_rotor_dont_mangle_8(genrand());
		
		// Test if norm2(V) is positive, otherwise do not perform the test.
		// (because with a negative norm2(v), the reverse is not the inverse)
		if (norm2_returns_scalar(V) <= 0.0) continue;

		tmp = unit(V);
		V = tmp;
		IV = versorInverse(V);

		// avoid near-singular versors
		if ((V.LargestCoordinate() > 2.0) ||
			(IV.LargestCoordinate() > 2.0))
			continue;
		 

		// generate random SMV 
		X = random_I3_t_dont_mangle_13(genrand());

		//  apply random versor to random SMV, convert to GMV
		VX = applyVersorWI(V, X, IV);
		gmvVX2 = new mv(VX);

		// convert all to GMV type, apply versor too as GMV
		gmvV = new mv(V);
//		gmvIV = new mv(IV);
		gmvX = new mv(X);
		gmvVX = applyVersor(gmvV, gmvX);
		
		// convert GMV back and forth to return type to fix possible constant coordinates
		tmpVX2 = new trivector(gmvVX);
		gmvVX = new mv(tmpVX2);
		
		// see if VX equals gmvVX
		dif = subtract(gmvVX, gmvVX2);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyVersorWI() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}


static int test_applyOM_dont_mangle_411(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	double[] OMmatrix = new double[3 * 3];
	int i, d, d2;
	int nbRandomVectors;
	om randomOM;
	mv[] randomVectors = new mv[3];
	mv[] OMrandomVectors = new mv[3];
	mv op1, op2, dif;
	int basisVectorBitmap = -1;
	int vectorGrade = 1;

	for (i = 0; i < NB_LOOPS; i++) {
		// init random outermorphism matrix
		for (d2 = 0; d2 < 3 * 3; d2++)
			OMmatrix[d2] = genrand();
		
		// init random OM
		randomOM = new om(OMmatrix);
		
		// get n < 3 random vectors stored in GMV
		nbRandomVectors = (int)(3.0 * genrand());
		if (nbRandomVectors < 1) nbRandomVectors = 1;
		for (d = 0; d < nbRandomVectors; d++) {
			randomVectors[d] = random_blade_dont_mangle_23_returns_mv(genrand(), vectorGrade, basisVectorBitmap);
			OMrandomVectors[d] = applyOM(randomOM, randomVectors[d]);
		}
		
		// compute outer product of randomVectors, OMrandomVectors
		op1 = randomVectors[0];
		op2 = OMrandomVectors[0];
		for (d = 1; d < nbRandomVectors; d++) {
			op1 = op(op1, randomVectors[d]);
			op2 = op(op2, OMrandomVectors[d]);
		}
		
		// apply OM to op1
		op1 = applyOM(randomOM, op1);
		
		// test equality
		dif = subtract(op1, op2);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("applyOM() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
		
	return 1; // success
}

static int test_applyOM_dont_mangle_412(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector[] rangeVectors = new vector[3];
	om refOM;
	om testOM;
	vector randomDomainSmv;
	vector rangeSmv, tmp;
	mv mv1, mv2, mv3, dif;
	int d, i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// init two range OMs (GOM and SOM)
		for (d = 0; d < 3; d++) {
			rangeVectors[d] = random_vector_dont_mangle_2(genrand());
		}
		refOM = new om(rangeVectors[0], rangeVectors[1], rangeVectors[2]);
		testOM = new om(refOM);
		
		
		// apply OM directly to domainSMV
		randomDomainSmv = random_vector_dont_mangle_2(genrand());
		
		rangeSmv = applyOM(testOM, randomDomainSmv);
		mv3 = new mv(rangeSmv);
		
		// convert domain SMV to GMV, apply to GMV
		mv1 = new mv(randomDomainSmv);
		mv2 = applyOM(refOM, mv1);
		
		// get rid of extra coordinates outside the range of the testOM
 
		tmp = new vector(mv2);
		mv2 = new mv(tmp);
 
		
		// test equality
		dif = subtract(mv2, mv3);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyOM() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
		
	return 1; // success
}

static int test_applyOM_dont_mangle_413(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector[] rangeVectors = new vector[3];
	om refOM;
	om testOM;
	bivector randomDomainSmv;
	bivector rangeSmv, tmp;
	mv mv1, mv2, mv3, dif;
	int d, i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// init two range OMs (GOM and SOM)
		for (d = 0; d < 3; d++) {
			rangeVectors[d] = random_vector_dont_mangle_2(genrand());
		}
		refOM = new om(rangeVectors[0], rangeVectors[1], rangeVectors[2]);
		testOM = new om(refOM);
		
		
		// apply OM directly to domainSMV
		randomDomainSmv = random_bivector_dont_mangle_4(genrand());
		
		rangeSmv = applyOM(testOM, randomDomainSmv);
		mv3 = new mv(rangeSmv);
		
		// convert domain SMV to GMV, apply to GMV
		mv1 = new mv(randomDomainSmv);
		mv2 = applyOM(refOM, mv1);
		
		// get rid of extra coordinates outside the range of the testOM
 
		tmp = new bivector(mv2);
		mv2 = new mv(tmp);
 
		
		// test equality
		dif = subtract(mv2, mv3);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyOM() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
		
	return 1; // success
}

static int test_applyOM_dont_mangle_414(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	vector[] rangeVectors = new vector[3];
	om refOM;
	om testOM;
	trivector randomDomainSmv;
	trivector rangeSmv, tmp;
	mv mv1, mv2, mv3, dif;
	int d, i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// init two range OMs (GOM and SOM)
		for (d = 0; d < 3; d++) {
			rangeVectors[d] = random_vector_dont_mangle_2(genrand());
		}
		refOM = new om(rangeVectors[0], rangeVectors[1], rangeVectors[2]);
		testOM = new om(refOM);
		
		
		// apply OM directly to domainSMV
		randomDomainSmv = random_trivector_dont_mangle_7(genrand());
		
		rangeSmv = applyOM(testOM, randomDomainSmv);
		mv3 = new mv(rangeSmv);
		
		// convert domain SMV to GMV, apply to GMV
		mv1 = new mv(randomDomainSmv);
		mv2 = applyOM(refOM, mv1);
		
		// get rid of extra coordinates outside the range of the testOM
 
		tmp = new trivector(mv2);
		mv2 = new mv(tmp);
 
		
		// test equality
		dif = subtract(mv2, mv3);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyOM() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
		
	return 1; // success
}

static int test_applyOM_dont_mangle_415(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector[] rangeVectors = new vector[3];
	om refOM;
	grade1OM testOM;
	vector randomDomainSmv;
	vector rangeSmv, tmp;
	mv mv1, mv2, mv3, dif;
	int d, i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// init two range OMs (GOM and SOM)
		for (d = 0; d < 3; d++) {
			rangeVectors[d] = random_vector_dont_mangle_2(genrand());
		}
		refOM = new om(rangeVectors[0], rangeVectors[1], rangeVectors[2]);
		testOM = new grade1OM(refOM);
		refOM = new om(testOM);
		
		
		// apply OM directly to domainSMV
		randomDomainSmv = random_vector_dont_mangle_2(genrand());
		
		rangeSmv = applyOM(testOM, randomDomainSmv);
		mv3 = new mv(rangeSmv);
		
		// convert domain SMV to GMV, apply to GMV
		mv1 = new mv(randomDomainSmv);
		mv2 = applyOM(refOM, mv1);
		
		// get rid of extra coordinates outside the range of the testOM
 
		tmp = new vector(mv2);
		mv2 = new mv(tmp);
 
		
		// test equality
		dif = subtract(mv2, mv3);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyOM() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
		
	return 1; // success
}

static int test_applyOM_dont_mangle_416(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector[] rangeVectors = new vector[3];
	om refOM;
	grade2OM testOM;
	bivector randomDomainSmv;
	bivector rangeSmv, tmp;
	mv mv1, mv2, mv3, dif;
	int d, i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// init two range OMs (GOM and SOM)
		for (d = 0; d < 3; d++) {
			rangeVectors[d] = random_vector_dont_mangle_2(genrand());
		}
		refOM = new om(rangeVectors[0], rangeVectors[1], rangeVectors[2]);
		testOM = new grade2OM(refOM);
		refOM = new om(testOM);
		
		
		// apply OM directly to domainSMV
		randomDomainSmv = random_bivector_dont_mangle_4(genrand());
		
		rangeSmv = applyOM(testOM, randomDomainSmv);
		mv3 = new mv(rangeSmv);
		
		// convert domain SMV to GMV, apply to GMV
		mv1 = new mv(randomDomainSmv);
		mv2 = applyOM(refOM, mv1);
		
		// get rid of extra coordinates outside the range of the testOM
 
		tmp = new bivector(mv2);
		mv2 = new mv(tmp);
 
		
		// test equality
		dif = subtract(mv2, mv3);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyOM() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
		
	return 1; // success
}

static int test_applyOM_dont_mangle_417(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	vector[] rangeVectors = new vector[3];
	om refOM;
	grade3OM testOM;
	trivector randomDomainSmv;
	trivector rangeSmv, tmp;
	mv mv1, mv2, mv3, dif;
	int d, i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// init two range OMs (GOM and SOM)
		for (d = 0; d < 3; d++) {
			rangeVectors[d] = random_vector_dont_mangle_2(genrand());
		}
		refOM = new om(rangeVectors[0], rangeVectors[1], rangeVectors[2]);
		testOM = new grade3OM(refOM);
		refOM = new om(testOM);
		
		
		// apply OM directly to domainSMV
		randomDomainSmv = random_trivector_dont_mangle_7(genrand());
		
		rangeSmv = applyOM(testOM, randomDomainSmv);
		mv3 = new mv(rangeSmv);
		
		// convert domain SMV to GMV, apply to GMV
		mv1 = new mv(randomDomainSmv);
		mv2 = applyOM(refOM, mv1);
		
		// get rid of extra coordinates outside the range of the testOM
 
		tmp = new trivector(mv2);
		mv2 = new mv(tmp);
 
		
		// test equality
		dif = subtract(mv2, mv3);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyOM() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
		
	return 1; // success
}

static int test_applyOM_dont_mangle_418(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	vector[] rangeVectors = new vector[3];
	om refOM;
	grade1OM testOM;
	e1_t randomDomainSmv;
	vector rangeSmv, tmp;
	mv mv1, mv2, mv3, dif;
	int d, i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// init two range OMs (GOM and SOM)
		for (d = 0; d < 3; d++) {
			rangeVectors[d] = random_vector_dont_mangle_2(genrand());
		}
		refOM = new om(rangeVectors[0], rangeVectors[1], rangeVectors[2]);
		testOM = new grade1OM(refOM);
		refOM = new om(testOM);
		
		
		// apply OM directly to domainSMV
		randomDomainSmv = random_e1_t_dont_mangle_10(genrand());
		
		rangeSmv = applyOM(testOM, randomDomainSmv);
		mv3 = new mv(rangeSmv);
		
		// convert domain SMV to GMV, apply to GMV
		mv1 = new mv(randomDomainSmv);
		mv2 = applyOM(refOM, mv1);
		
		// get rid of extra coordinates outside the range of the testOM
 
		tmp = new vector(mv2);
		mv2 = new mv(tmp);
 
		
		// test equality
		dif = subtract(mv2, mv3);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyOM() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
		
	return 1; // success
}

static int test_applyOM_dont_mangle_419(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	vector[] rangeVectors = new vector[3];
	om refOM;
	grade1OM testOM;
	e2_t randomDomainSmv;
	vector rangeSmv, tmp;
	mv mv1, mv2, mv3, dif;
	int d, i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// init two range OMs (GOM and SOM)
		for (d = 0; d < 3; d++) {
			rangeVectors[d] = random_vector_dont_mangle_2(genrand());
		}
		refOM = new om(rangeVectors[0], rangeVectors[1], rangeVectors[2]);
		testOM = new grade1OM(refOM);
		refOM = new om(testOM);
		
		
		// apply OM directly to domainSMV
		randomDomainSmv = random_e2_t_dont_mangle_11(genrand());
		
		rangeSmv = applyOM(testOM, randomDomainSmv);
		mv3 = new mv(rangeSmv);
		
		// convert domain SMV to GMV, apply to GMV
		mv1 = new mv(randomDomainSmv);
		mv2 = applyOM(refOM, mv1);
		
		// get rid of extra coordinates outside the range of the testOM
 
		tmp = new vector(mv2);
		mv2 = new mv(tmp);
 
		
		// test equality
		dif = subtract(mv2, mv3);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyOM() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
		
	return 1; // success
}

static int test_applyOM_dont_mangle_420(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	vector[] rangeVectors = new vector[3];
	om refOM;
	grade1OM testOM;
	e3_t randomDomainSmv;
	vector rangeSmv, tmp;
	mv mv1, mv2, mv3, dif;
	int d, i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// init two range OMs (GOM and SOM)
		for (d = 0; d < 3; d++) {
			rangeVectors[d] = random_vector_dont_mangle_2(genrand());
		}
		refOM = new om(rangeVectors[0], rangeVectors[1], rangeVectors[2]);
		testOM = new grade1OM(refOM);
		refOM = new om(testOM);
		
		
		// apply OM directly to domainSMV
		randomDomainSmv = random_e3_t_dont_mangle_75(genrand());
		
		rangeSmv = applyOM(testOM, randomDomainSmv);
		mv3 = new mv(rangeSmv);
		
		// convert domain SMV to GMV, apply to GMV
		mv1 = new mv(randomDomainSmv);
		mv2 = applyOM(refOM, mv1);
		
		// get rid of extra coordinates outside the range of the testOM
 
		tmp = new vector(mv2);
		mv2 = new mv(tmp);
 
		
		// test equality
		dif = subtract(mv2, mv3);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyOM() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
		
	return 1; // success
}

static int test_applyOM_dont_mangle_421(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	vector[] rangeVectors = new vector[3];
	om refOM;
	grade3OM testOM;
	I3_t randomDomainSmv;
	trivector rangeSmv, tmp;
	mv mv1, mv2, mv3, dif;
	int d, i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// init two range OMs (GOM and SOM)
		for (d = 0; d < 3; d++) {
			rangeVectors[d] = random_vector_dont_mangle_2(genrand());
		}
		refOM = new om(rangeVectors[0], rangeVectors[1], rangeVectors[2]);
		testOM = new grade3OM(refOM);
		refOM = new om(testOM);
		
		
		// apply OM directly to domainSMV
		randomDomainSmv = random_I3_t_dont_mangle_13(genrand());
		
		rangeSmv = applyOM(testOM, randomDomainSmv);
		mv3 = new mv(rangeSmv);
		
		// convert domain SMV to GMV, apply to GMV
		mv1 = new mv(randomDomainSmv);
		mv2 = applyOM(refOM, mv1);
		
		// get rid of extra coordinates outside the range of the testOM
 
		tmp = new trivector(mv2);
		mv2 = new mv(tmp);
 
		
		// test equality
		dif = subtract(mv2, mv3);
		if (dif.LargestCoordinate() > (1E-12 )) {
			Console.WriteLine("applyOM() test failed (largest coordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
		
	return 1; // success
}

static int test_div_dont_mangle_422(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	int i;
	mv A, B, C, dif;
	double divider;
	
	int basisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_blade_dont_mangle_23_returns_mv(genrand(), (int)(genrand() * 3.5), basisVectorBitmap);
		
		divider = 0.01 + genrand();
		
		B = div(A, divider);
		
		C = gp(B, divider);
		
		// see if C equals A
		dif = subtract(C, A);
		if (dif.LargestCoordinate() > (1E-14 )) {
			Console.WriteLine("div() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}

static int test_div_dont_mangle_423(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	int i;
	vector A;
	vector B;
	mv gmvA, gmvB, C, dif;
	double divider;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random smv
		A = random_vector_dont_mangle_2(genrand());
		
		divider = 0.01 + genrand();
		
		B = div(A, divider);
		
		gmvB = new mv(B);
		C = gp(gmvB, divider);
		
		gmvA = new mv(A);
		
		// see if C equals A
		dif = subtract(C, gmvA);
		if (dif.LargestCoordinate() > (1E-13 )) {
			Console.WriteLine("div() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}

static int test_div_dont_mangle_424(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	int i;
	bivector A;
	bivector B;
	mv gmvA, gmvB, C, dif;
	double divider;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random smv
		A = random_bivector_dont_mangle_4(genrand());
		
		divider = 0.01 + genrand();
		
		B = div(A, divider);
		
		gmvB = new mv(B);
		C = gp(gmvB, divider);
		
		gmvA = new mv(A);
		
		// see if C equals A
		dif = subtract(C, gmvA);
		if (dif.LargestCoordinate() > (1E-13 )) {
			Console.WriteLine("div() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}

static int test_div_dont_mangle_425(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	int i;
	trivector A;
	trivector B;
	mv gmvA, gmvB, C, dif;
	double divider;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random smv
		A = random_trivector_dont_mangle_7(genrand());
		
		divider = 0.01 + genrand();
		
		B = div(A, divider);
		
		gmvB = new mv(B);
		C = gp(gmvB, divider);
		
		gmvA = new mv(A);
		
		// see if C equals A
		dif = subtract(C, gmvA);
		if (dif.LargestCoordinate() > (1E-13 )) {
			Console.WriteLine("div() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}

static int test_div_dont_mangle_426(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	int i;
	rotor A;
	rotor B;
	mv gmvA, gmvB, C, dif;
	double divider;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random smv
		A = random_rotor_dont_mangle_8(genrand());
		
		divider = 0.01 + genrand();
		
		B = div(A, divider);
		
		gmvB = new mv(B);
		C = gp(gmvB, divider);
		
		gmvA = new mv(A);
		
		// see if C equals A
		dif = subtract(C, gmvA);
		if (dif.LargestCoordinate() > (1E-13 )) {
			Console.WriteLine("div() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}

static int test_div_dont_mangle_427(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	int i;
	e1_t A;
	vector B;
	mv gmvA, gmvB, C, dif;
	double divider;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random smv
		A = random_e1_t_dont_mangle_10(genrand());
		
		divider = 0.01 + genrand();
		
		B = div(A, divider);
		
		gmvB = new mv(B);
		C = gp(gmvB, divider);
		
		gmvA = new mv(A);
		
		// see if C equals A
		dif = subtract(C, gmvA);
		if (dif.LargestCoordinate() > (1E-13 )) {
			Console.WriteLine("div() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}

static int test_div_dont_mangle_428(int NB_TESTS_SCALER) {
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	int i;
	I3_t A;
	trivector B;
	mv gmvA, gmvB, C, dif;
	double divider;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random smv
		A = random_I3_t_dont_mangle_13(genrand());
		
		divider = 0.01 + genrand();
		
		B = div(A, divider);
		
		gmvB = new mv(B);
		C = gp(gmvB, divider);
		
		gmvA = new mv(A);
		
		// see if C equals A
		dif = subtract(C, gmvA);
		if (dif.LargestCoordinate() > (1E-13 )) {
			Console.WriteLine("div() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}

	return 1; // success
}

static int test_dual_dont_mangle_429(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, C, dif;
	int i;
	int basisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_blade_dont_mangle_23_returns_mv(genrand(), (int)(genrand() * 3.5), basisVectorBitmap);
		
		B = dual(A);
		
		C = undual(B);
		
		// check if equal to original:
		dif = subtract(A, C);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("dual() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_undual_dont_mangle_430(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, C, dif;
	int i;
	int basisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_blade_dont_mangle_23_returns_mv(genrand(), (int)(genrand() * 3.5), basisVectorBitmap);
		
		B = undual(A);
		
		C = dual(B);
		
		// check if equal to original:
		dif = subtract(A, C);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("undual() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_dual_dont_mangle_433(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	bivector B;
	mv gmvA, gmvB, C, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_vector_dont_mangle_2(genrand());
		gmvA = new mv(A);
		
		B = dual(A);
		gmvB = new mv(B);
		
		C = undual(gmvB);
		
		// check if equal to original:
		dif = subtract(gmvA, C);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("dual() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_undual_dont_mangle_434(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	bivector B;
	mv gmvA, gmvB, C, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_vector_dont_mangle_2(genrand());
		gmvA = new mv(A);
		
		B = undual(A);
		gmvB = new mv(B);
		
		C = dual(gmvB);
		
		// check if equal to original:
		dif = subtract(gmvA, C);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("undual() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_dual_dont_mangle_435(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	vector B;
	mv gmvA, gmvB, C, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_bivector_dont_mangle_4(genrand());
		gmvA = new mv(A);
		
		B = dual(A);
		gmvB = new mv(B);
		
		C = undual(gmvB);
		
		// check if equal to original:
		dif = subtract(gmvA, C);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("dual() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_undual_dont_mangle_436(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	vector B;
	mv gmvA, gmvB, C, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_bivector_dont_mangle_4(genrand());
		gmvA = new mv(A);
		
		B = undual(A);
		gmvB = new mv(B);
		
		C = dual(gmvB);
		
		// check if equal to original:
		dif = subtract(gmvA, C);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("undual() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_dual_dont_mangle_437(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	oddVersor B;
	mv gmvA, gmvB, C, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_rotor_dont_mangle_8(genrand());
		gmvA = new mv(A);
		
		B = dual(A);
		gmvB = new mv(B);
		
		C = undual(gmvB);
		
		// check if equal to original:
		dif = subtract(gmvA, C);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("dual() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_undual_dont_mangle_438(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	oddVersor B;
	mv gmvA, gmvB, C, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_rotor_dont_mangle_8(genrand());
		gmvA = new mv(A);
		
		B = undual(A);
		gmvB = new mv(B);
		
		C = dual(gmvB);
		
		// check if equal to original:
		dif = subtract(gmvA, C);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("undual() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_dual_dont_mangle_439(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	oddVersor A;
	rotor B;
	mv gmvA, gmvB, C, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_oddVersor_dont_mangle_93(genrand());
		gmvA = new mv(A);
		
		B = dual(A);
		gmvB = new mv(B);
		
		C = undual(gmvB);
		
		// check if equal to original:
		dif = subtract(gmvA, C);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("dual() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_undual_dont_mangle_440(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	oddVersor A;
	rotor B;
	mv gmvA, gmvB, C, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_oddVersor_dont_mangle_93(genrand());
		gmvA = new mv(A);
		
		B = undual(A);
		gmvB = new mv(B);
		
		C = dual(gmvB);
		
		// check if equal to original:
		dif = subtract(gmvA, C);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("undual() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_dual_dont_mangle_441(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	trivector A;
	double B;
	mv gmvA, gmvB, C, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_trivector_dont_mangle_7(genrand());
		gmvA = new mv(A);
		
		B = dual(A);
		gmvB = new mv(B);
		
		C = undual(gmvB);
		
		// check if equal to original:
		dif = subtract(gmvA, C);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("dual() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_undual_dont_mangle_442(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	trivector A;
	double B;
	mv gmvA, gmvB, C, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_trivector_dont_mangle_7(genrand());
		gmvA = new mv(A);
		
		B = undual(A);
		gmvB = new mv(B);
		
		C = dual(gmvB);
		
		// check if equal to original:
		dif = subtract(gmvA, C);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("undual() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_dual_dont_mangle_443(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	e1_t A;
	bivector B;
	mv gmvA, gmvB, C, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_e1_t_dont_mangle_10(genrand());
		gmvA = new mv(A);
		
		B = dual(A);
		gmvB = new mv(B);
		
		C = undual(gmvB);
		
		// check if equal to original:
		dif = subtract(gmvA, C);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("dual() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_undual_dont_mangle_444(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	e2_t A;
	bivector B;
	mv gmvA, gmvB, C, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_e2_t_dont_mangle_11(genrand());
		gmvA = new mv(A);
		
		B = undual(A);
		gmvB = new mv(B);
		
		C = dual(gmvB);
		
		// check if equal to original:
		dif = subtract(gmvA, C);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("undual() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_dual_dont_mangle_445(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	I3_t A;
	double B;
	mv gmvA, gmvB, C, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_I3_t_dont_mangle_13(genrand());
		gmvA = new mv(A);
		
		B = dual(A);
		gmvB = new mv(B);
		
		C = undual(gmvB);
		
		// check if equal to original:
		dif = subtract(gmvA, C);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("dual() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_undual_dont_mangle_446(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	I3_t A;
	double B;
	mv gmvA, gmvB, C, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_I3_t_dont_mangle_13(genrand());
		gmvA = new mv(A);
		
		B = undual(A);
		gmvB = new mv(B);
		
		C = dual(gmvB);
		
		// check if equal to original:
		dif = subtract(gmvA, C);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("undual() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_equals_dont_mangle_447(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, C;
	double s, eps = 0.2;
	int g, i;
	bool e, l;
	int basisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		g = (int)(genrand() * 3.5);
		A = random_versor_dont_mangle_1_returns_mv(s, g, basisVectorBitmap);
		B = random_versor_dont_mangle_1_returns_mv(1.1 * eps, g, basisVectorBitmap);
		C = add(B, A);
		
		// check if equals thinks A if is equal to itself
		e = equals(A, A, eps);
		if (!e) {
			Console.WriteLine("equals() failed (variable not equal to itself)\n");
			return 0; // failure
		}
		
		// check if equals thinks A if is equal C
		e = equals(A, C, eps);
		
		// use mv_largestCoordinate() to verify
		l = !(B.LargestCoordinate() > eps);
		
		if (e != l) {
			Console.WriteLine("equals() test failed");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_equals_dont_mangle_448(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 6;
	vector A;
	vector B;
	mv gA, gB, gC;
	double s, eps = 0.5; // note the really large epsilon
	int i;
	bool e, l;
	bool e1, e2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_vector_dont_mangle_2(s);
		B = random_vector_dont_mangle_2(s);
		
		// check if equals thinks A if is equal to itself
		e1 = equals(A, A, eps);
		e2 = equals(B, B, eps);
		if ((!e1) || (!e2)) {
			Console.WriteLine("equals() failed (variable not equal to itself)\n");
			return 0; // failure
		}		// convert A and B to gmv
		gA = new mv(A);
		gB = new mv(B);
		
		// gC = gB - gA, get largest coord
		gC = subtract(gB, gA);
		// use largestCoordinate() to verify
		l = !(gC.LargestCoordinate() > eps);
		
		// check if equals thinks A if is equal B
		e = equals(A, B, eps);
		
		if (e != l) {
			Console.WriteLine("equals() test failed\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_equals_dont_mangle_449(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 6;
	bivector A;
	bivector B;
	mv gA, gB, gC;
	double s, eps = 0.5; // note the really large epsilon
	int i;
	bool e, l;
	bool e1, e2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_bivector_dont_mangle_4(s);
		B = random_bivector_dont_mangle_4(s);
		
		// check if equals thinks A if is equal to itself
		e1 = equals(A, A, eps);
		e2 = equals(B, B, eps);
		if ((!e1) || (!e2)) {
			Console.WriteLine("equals() failed (variable not equal to itself)\n");
			return 0; // failure
		}		// convert A and B to gmv
		gA = new mv(A);
		gB = new mv(B);
		
		// gC = gB - gA, get largest coord
		gC = subtract(gB, gA);
		// use largestCoordinate() to verify
		l = !(gC.LargestCoordinate() > eps);
		
		// check if equals thinks A if is equal B
		e = equals(A, B, eps);
		
		if (e != l) {
			Console.WriteLine("equals() test failed\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_equals_dont_mangle_450(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	rotor A;
	rotor B;
	mv gA, gB, gC;
	double s, eps = 0.5; // note the really large epsilon
	int i;
	bool e, l;
	bool e1, e2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_rotor_dont_mangle_8(s);
		B = random_rotor_dont_mangle_8(s);
		
		// check if equals thinks A if is equal to itself
		e1 = equals(A, A, eps);
		e2 = equals(B, B, eps);
		if ((!e1) || (!e2)) {
			Console.WriteLine("equals() failed (variable not equal to itself)\n");
			return 0; // failure
		}		// convert A and B to gmv
		gA = new mv(A);
		gB = new mv(B);
		
		// gC = gB - gA, get largest coord
		gC = subtract(gB, gA);
		// use largestCoordinate() to verify
		l = !(gC.LargestCoordinate() > eps);
		
		// check if equals thinks A if is equal B
		e = equals(A, B, eps);
		
		if (e != l) {
			Console.WriteLine("equals() test failed\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_equals_dont_mangle_451(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 7;
	bivector A;
	rotor B;
	mv gA, gB, gC;
	double s, eps = 0.5; // note the really large epsilon
	int i;
	bool e, l;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_bivector_dont_mangle_4(s);
		B = random_rotor_dont_mangle_8(s);
		
		// convert A and B to gmv
		gA = new mv(A);
		gB = new mv(B);
		
		// gC = gB - gA, get largest coord
		gC = subtract(gB, gA);
		// use largestCoordinate() to verify
		l = !(gC.LargestCoordinate() > eps);
		
		// check if equals thinks A if is equal B
		e = equals(A, B, eps);
		
		if (e != l) {
			Console.WriteLine("equals() test failed\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_equals_dont_mangle_452(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 2;
	trivector A;
	trivector B;
	mv gA, gB, gC;
	double s, eps = 0.5; // note the really large epsilon
	int i;
	bool e, l;
	bool e1, e2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_trivector_dont_mangle_7(s);
		B = random_trivector_dont_mangle_7(s);
		
		// check if equals thinks A if is equal to itself
		e1 = equals(A, A, eps);
		e2 = equals(B, B, eps);
		if ((!e1) || (!e2)) {
			Console.WriteLine("equals() failed (variable not equal to itself)\n");
			return 0; // failure
		}		// convert A and B to gmv
		gA = new mv(A);
		gB = new mv(B);
		
		// gC = gB - gA, get largest coord
		gC = subtract(gB, gA);
		// use largestCoordinate() to verify
		l = !(gC.LargestCoordinate() > eps);
		
		// check if equals thinks A if is equal B
		e = equals(A, B, eps);
		
		if (e != l) {
			Console.WriteLine("equals() test failed\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_equals_dont_mangle_453(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 7;
	rotor A;
	bivector B;
	mv gA, gB, gC;
	double s, eps = 0.5; // note the really large epsilon
	int i;
	bool e, l;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_rotor_dont_mangle_8(s);
		B = random_bivector_dont_mangle_4(s);
		
		// convert A and B to gmv
		gA = new mv(A);
		gB = new mv(B);
		
		// gC = gB - gA, get largest coord
		gC = subtract(gB, gA);
		// use largestCoordinate() to verify
		l = !(gC.LargestCoordinate() > eps);
		
		// check if equals thinks A if is equal B
		e = equals(A, B, eps);
		
		if (e != l) {
			Console.WriteLine("equals() test failed\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_equals_dont_mangle_454(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 2;
	e1_t A;
	e1_t B;
	mv gA, gB, gC;
	double s, eps = 0.5; // note the really large epsilon
	int i;
	bool e, l;
	bool e1, e2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_e1_t_dont_mangle_10(s);
		B = random_e1_t_dont_mangle_10(s);
		
		// check if equals thinks A if is equal to itself
		e1 = equals(A, A, eps);
		e2 = equals(B, B, eps);
		if ((!e1) || (!e2)) {
			Console.WriteLine("equals() failed (variable not equal to itself)\n");
			return 0; // failure
		}		// convert A and B to gmv
		gA = new mv(A);
		gB = new mv(B);
		
		// gC = gB - gA, get largest coord
		gC = subtract(gB, gA);
		// use largestCoordinate() to verify
		l = !(gC.LargestCoordinate() > eps);
		
		// check if equals thinks A if is equal B
		e = equals(A, B, eps);
		
		if (e != l) {
			Console.WriteLine("equals() test failed\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_equals_dont_mangle_455(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 2;
	e2_t A;
	I3_t B;
	mv gA, gB, gC;
	double s, eps = 0.5; // note the really large epsilon
	int i;
	bool e, l;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_e2_t_dont_mangle_11(s);
		B = random_I3_t_dont_mangle_13(s);
		
		// convert A and B to gmv
		gA = new mv(A);
		gB = new mv(B);
		
		// gC = gB - gA, get largest coord
		gC = subtract(gB, gA);
		// use largestCoordinate() to verify
		l = !(gC.LargestCoordinate() > eps);
		
		// check if equals thinks A if is equal B
		e = equals(A, B, eps);
		
		if (e != l) {
			Console.WriteLine("equals() test failed\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_extractGrade_dont_mangle_456(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, C, G0, G1, G2, G3;
	int i;
	int basisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		A = random_versor_dont_mangle_1_returns_mv(genrand(), (int)(genrand() * 3.5), basisVectorBitmap);
		// split it up into groups/grades:
		G0 = extractGrade(A, GroupBitmap.GROUP_0);
		G1 = extractGrade(A, GroupBitmap.GROUP_1);
		G2 = extractGrade(A, GroupBitmap.GROUP_2);
		G3 = extractGrade(A, GroupBitmap.GROUP_3);
		// sum all into 'B'
		B = new mv();
		B = add(B, G0);
		B = add(B, G1);
		B = add(B, G2);
		B = add(B, G3);

		// check if equal to original:
		C = subtract(A, B);
		if (C.LargestCoordinate() > 1E-14) {
			Console.WriteLine("extractGrade() test failed (largestCoordinate = " + (double)C.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_extractGrade2_dont_mangle_457(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, C, D;
	int i;
	int basisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		A = random_versor_dont_mangle_1_returns_mv(genrand(), (int)(genrand() * 3.5), basisVectorBitmap);
		
		B = extractGrade2(A);
		
		C = extractGrade(A, 0 | GroupBitmap.GROUP_0 | GroupBitmap.GROUP_1 | GroupBitmap.GROUP_3);
		
		// sum all into 'B'
		D = add(B, C);

		// check if equal to original:
		C = subtract(A, D);
		if (C.LargestCoordinate() > 1E-14) {
			Console.WriteLine("extractGrade2() test failed");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_extractGrade0_dont_mangle_458(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 5;
	rotor A;
	mv gA, gB, gC, gD;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		A = random_rotor_dont_mangle_8(genrand());
		
		gB = new mv(extractGrade0(A));

		gA = new mv(A);
		gC = extractGrade(gA, 0 | GroupBitmap.GROUP_0);
		
		// check if equal to original:
		gD = subtract(gB, gC);
		if (gD.LargestCoordinate() > 1E-14) {
			Console.WriteLine("extractGrade0() test failed");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_extractGrade2_dont_mangle_459(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 5;
	rotor A;
	mv gA, gB, gC, gD;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		A = random_rotor_dont_mangle_8(genrand());
		
		gB = new mv(extractGrade2(A));

		gA = new mv(A);
		gC = extractGrade(gA, 0 | GroupBitmap.GROUP_2);
		
		// check if equal to original:
		gD = subtract(gB, gC);
		if (gD.LargestCoordinate() > 1E-14) {
			Console.WriteLine("extractGrade2() test failed");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_extractGrade0_dont_mangle_460(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 5;
	oddVersor A;
	mv gA, gB, gC, gD;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		A = random_oddVersor_dont_mangle_93(genrand());
		
		gB = new mv(extractGrade0(A));

		gA = new mv(A);
		gC = extractGrade(gA, 0 | GroupBitmap.GROUP_0);
		
		// check if equal to original:
		gD = subtract(gB, gC);
		if (gD.LargestCoordinate() > 1E-14) {
			Console.WriteLine("extractGrade0() test failed");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_extractGrade1_dont_mangle_461(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 5;
	oddVersor A;
	mv gA, gB, gC, gD;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		A = random_oddVersor_dont_mangle_93(genrand());
		
		gB = new mv(extractGrade1(A));

		gA = new mv(A);
		gC = extractGrade(gA, 0 | GroupBitmap.GROUP_1);
		
		// check if equal to original:
		gD = subtract(gB, gC);
		if (gD.LargestCoordinate() > 1E-14) {
			Console.WriteLine("extractGrade1() test failed");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_extractGrade2_dont_mangle_462(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 5;
	oddVersor A;
	mv gA, gB, gC, gD;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		A = random_oddVersor_dont_mangle_93(genrand());
		
		gB = new mv(extractGrade2(A));

		gA = new mv(A);
		gC = extractGrade(gA, 0 | GroupBitmap.GROUP_2);
		
		// check if equal to original:
		gD = subtract(gB, gC);
		if (gD.LargestCoordinate() > 1E-14) {
			Console.WriteLine("extractGrade2() test failed");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_extractGrade3_dont_mangle_463(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 5;
	oddVersor A;
	mv gA, gB, gC, gD;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		A = random_oddVersor_dont_mangle_93(genrand());
		
		gB = new mv(extractGrade3(A));

		gA = new mv(A);
		gC = extractGrade(gA, 0 | GroupBitmap.GROUP_3);
		
		// check if equal to original:
		gD = subtract(gB, gC);
		if (gD.LargestCoordinate() > 1E-14) {
			Console.WriteLine("extractGrade3() test failed");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_extractGrade0_dont_mangle_464(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 2;
	e1_t A;
	mv gA, gB, gC, gD;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		A = random_e1_t_dont_mangle_10(genrand());
		
		gB = new mv(extractGrade0(A));

		gA = new mv(A);
		gC = extractGrade(gA, 0 | GroupBitmap.GROUP_0);
		
		// check if equal to original:
		gD = subtract(gB, gC);
		if (gD.LargestCoordinate() > 1E-14) {
			Console.WriteLine("extractGrade0() test failed");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_extractGrade1_dont_mangle_465(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 2;
	e2_t A;
	mv gA, gB, gC, gD;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		A = random_e2_t_dont_mangle_11(genrand());
		
		gB = new mv(extractGrade1(A));

		gA = new mv(A);
		gC = extractGrade(gA, 0 | GroupBitmap.GROUP_1);
		
		// check if equal to original:
		gD = subtract(gB, gC);
		if (gD.LargestCoordinate() > 1E-14) {
			Console.WriteLine("extractGrade1() test failed");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_extractGrade2_dont_mangle_466(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 2;
	e3_t A;
	mv gA, gB, gC, gD;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		A = random_e3_t_dont_mangle_75(genrand());
		
		gB = new mv(extractGrade2(A));

		gA = new mv(A);
		gC = extractGrade(gA, 0 | GroupBitmap.GROUP_2);
		
		// check if equal to original:
		gD = subtract(gB, gC);
		if (gD.LargestCoordinate() > 1E-14) {
			Console.WriteLine("extractGrade2() test failed");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_extractGrade3_dont_mangle_467(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 2;
	e1_t A;
	mv gA, gB, gC, gD;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		A = random_e1_t_dont_mangle_10(genrand());
		
		gB = new mv(extractGrade3(A));

		gA = new mv(A);
		gC = extractGrade(gA, 0 | GroupBitmap.GROUP_3);
		
		// check if equal to original:
		gD = subtract(gB, gC);
		if (gD.LargestCoordinate() > 1E-14) {
			Console.WriteLine("extractGrade3() test failed");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_extractGrade0_dont_mangle_468(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 2;
	I3_t A;
	mv gA, gB, gC, gD;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		A = random_I3_t_dont_mangle_13(genrand());
		
		gB = new mv(extractGrade0(A));

		gA = new mv(A);
		gC = extractGrade(gA, 0 | GroupBitmap.GROUP_0);
		
		// check if equal to original:
		gD = subtract(gB, gC);
		if (gD.LargestCoordinate() > 1E-14) {
			Console.WriteLine("extractGrade0() test failed");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_extractGrade1_dont_mangle_469(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 2;
	I3_t A;
	mv gA, gB, gC, gD;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		A = random_I3_t_dont_mangle_13(genrand());
		
		gB = new mv(extractGrade1(A));

		gA = new mv(A);
		gC = extractGrade(gA, 0 | GroupBitmap.GROUP_1);
		
		// check if equal to original:
		gD = subtract(gB, gC);
		if (gD.LargestCoordinate() > 1E-14) {
			Console.WriteLine("extractGrade1() test failed");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_extractGrade2_dont_mangle_470(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 2;
	I3_t A;
	mv gA, gB, gC, gD;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		A = random_I3_t_dont_mangle_13(genrand());
		
		gB = new mv(extractGrade2(A));

		gA = new mv(A);
		gC = extractGrade(gA, 0 | GroupBitmap.GROUP_2);
		
		// check if equal to original:
		gD = subtract(gB, gC);
		if (gD.LargestCoordinate() > 1E-14) {
			Console.WriteLine("extractGrade2() test failed");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_extractGrade3_dont_mangle_471(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 2;
	I3_t A;
	mv gA, gB, gC, gD;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		A = random_I3_t_dont_mangle_13(genrand());
		
		gB = new mv(extractGrade3(A));

		gA = new mv(A);
		gC = extractGrade(gA, 0 | GroupBitmap.GROUP_3);
		
		// check if equal to original:
		gD = subtract(gB, gC);
		if (gD.LargestCoordinate() > 1E-14) {
			Console.WriteLine("extractGrade3() test failed");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_gp_dont_mangle_472(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 32;
	mv A, B, C, D, E, V1, V2;
	int i;
	int o;
	int basisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		o = (genrand() < 0.5) ? 0 : 1; // even or odd?
		A = random_versor_dont_mangle_1_returns_mv(genrand(), ((int)(genrand() * 3.5) & 0xFFFE) + o, basisVectorBitmap);
		B = random_versor_dont_mangle_1_returns_mv(genrand(), ((int)(genrand() * 3.5) & 0xFFFE) + o, basisVectorBitmap);
		C = random_versor_dont_mangle_1_returns_mv(genrand(), ((int)(genrand() * 3.5) & 0xFFFE) + o, basisVectorBitmap);
		
		{ // test (A+B) C = A C + B C
			// D = A + B
			D = add(A, B);
			// V1 = D C
			V1 = gp(D, C);
			// D = A C
			D = gp(A, C);
			// E = B C
			E = gp(B, C);
			// V2 = D + E
			V2 = add(D, E);
			// test equality
			D = subtract(V1, V2);
			// use mv_largestCoordinate() to verify
			if (D.LargestCoordinate() > 1E-11) {
				Console.WriteLine("gp() test failed on '(A+B) C = A C + B C' (dif=" + D.LargestCoordinate() + ")");
				return 0; // failure
			}
		}
		
		{ // test A (B+C) = A B + A C
			// D = B + C
			D = add(B, C);
			// V1 = A D
			V1 = gp(A, D);
			// D = A B
			D = gp(A, B);
			// E = A C
			E = gp(A, C);
			// V2 = D + E
			V2 = add(D, E);
			// test equality
			D = subtract(V1, V2);
			// use largestCoordinate() to verify
			if (D.LargestCoordinate() > 1E-12) {
				Console.WriteLine("gp() test failed on 'A (B+C) = A B + A C' (dif=" + D.LargestCoordinate() + ")");
				return 0; // failure
			}
		}
		
		{ // test A (B C) = (A B) C
			// D = B C
			D = gp(B, C);
			// V1 = A D
			V1 = gp(A, D);
			// D = A B
			D = gp(A, B);
			// V2 = D C
			V2 = gp(D, C);
			// test equality
			D = subtract(V1, V2);
			// use largestCoordinate() to verify
			if (D.LargestCoordinate() > 1E-12) {
				Console.WriteLine("gp() test failed on 'A (B C) = (A B) C' (dif=" + D.LargestCoordinate() + ")");
				return 0; // failure
			}
		}		
	}
	return 1; // success
}

static int test_gp_dont_mangle_477(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 6;
	vector A;
	vector B;
	rotor C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_vector_dont_mangle_2(s);
		B = random_vector_dont_mangle_2(s);
		
		// A gp B
		C = gp(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = gp(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("gp() test failed (largestCoordinate = " + gA.LargestCoordinate() + ")");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_gp_dont_mangle_475(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 7;
	rotor A;
	vector B;
	oddVersor C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_rotor_dont_mangle_8(s);
		B = random_vector_dont_mangle_2(s);
		
		// A gp B
		C = gp(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = gp(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("gp() test failed (largestCoordinate = " + gA.LargestCoordinate() + ")");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_gp_dont_mangle_478(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 7;
	vector A;
	rotor B;
	oddVersor C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_vector_dont_mangle_2(s);
		B = random_rotor_dont_mangle_8(s);
		
		// A gp B
		C = gp(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = gp(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("gp() test failed (largestCoordinate = " + gA.LargestCoordinate() + ")");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_gp_dont_mangle_479(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	rotor A;
	rotor B;
	rotor C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_rotor_dont_mangle_8(s);
		B = random_rotor_dont_mangle_8(s);
		
		// A gp B
		C = gp(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = gp(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("gp() test failed (largestCoordinate = " + gA.LargestCoordinate() + ")");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_gp_dont_mangle_474(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 6;
	bivector A;
	bivector B;
	rotor C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_bivector_dont_mangle_4(s);
		B = random_bivector_dont_mangle_4(s);
		
		// A gp B
		C = gp(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = gp(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("gp() test failed (largestCoordinate = " + gA.LargestCoordinate() + ")");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_gp_dont_mangle_476(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 5;
	e1_t A;
	rotor B;
	oddVersor C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_e1_t_dont_mangle_10(s);
		B = random_rotor_dont_mangle_8(s);
		
		// A gp B
		C = gp(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = gp(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("gp() test failed (largestCoordinate = " + gA.LargestCoordinate() + ")");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_gp_dont_mangle_473(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 5;
	I3_t A;
	rotor B;
	oddVersor C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_I3_t_dont_mangle_13(s);
		B = random_rotor_dont_mangle_8(s);
		
		// A gp B
		C = gp(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = gp(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("gp() test failed (largestCoordinate = " + gA.LargestCoordinate() + ")");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_gp_dont_mangle_480(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	bivector A;
	e1_t B;
	oddVersor C;
	mv gA, gB, gC1, gC2;
	
	double s;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_bivector_dont_mangle_4(s);
		B = random_e1_t_dont_mangle_10(s);
		
		// A gp B
		C = gp(A, B);
		gC1 = new mv(C);
		
		// convert all A and B to gmv and add/subtract as GMVs
		gA = new mv(A);
		gB = new mv(B);
		gC2 = gp(gA, gB);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-13) {
			Console.WriteLine("gp() test failed (largestCoordinate = " + gA.LargestCoordinate() + ")");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_gradeBitmap_dont_mangle_481(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, tmp, randomBlade;
	int i, j;
	int basisVectorBitmap = -1;
	int gradeBitmap1, gradeBitmap2, nbBlades, grade;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get sum of random blades, keep track of grades used
		gradeBitmap1 = 0;
		A = new mv();
		nbBlades = (int)(genrand() * 3.5);
		for (j = 0; j < nbBlades; j++) {
			grade = (int)(genrand() * 3.5);
			gradeBitmap1 |= 1 << grade;
			randomBlade = random_blade_dont_mangle_23_returns_mv(1.0, grade, basisVectorBitmap);
			tmp = add(A, randomBlade);
			A = tmp;
		}
		
		gradeBitmap2 = gradeBitmap(A, 0.0);
		
		// check if grade bitmaps match
		if (gradeBitmap1 != gradeBitmap2) {
			Console.WriteLine("gradeBitmap() test failed (grade bitmap " + gradeBitmap1 + " vs " + gradeBitmap2 + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_gradeBitmap_dont_mangle_482(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	mv gmvA;
	rotor A;
	int i;
	int gradeBitmap1, gradeBitmap2;
	double threshold;
	
	for (i = 0; i < NB_LOOPS; i++) {
		threshold = 0.01 * genrand();
		A = random_rotor_dont_mangle_8(1.0);
		
		gradeBitmap1 = gradeBitmap(A, threshold);
		
		gmvA = new mv(A);
		
		gradeBitmap2 = gradeBitmap(gmvA, threshold);
		
		// check if grade bitmaps match
		if (gradeBitmap1 != gradeBitmap2) {
			Console.WriteLine("gradeBitmap() test failed (grade bitmap " + gradeBitmap1 + " vs " + gradeBitmap2 + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_gradeBitmap_dont_mangle_483(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	mv gmvA;
	vector A;
	int i;
	int gradeBitmap1, gradeBitmap2;
	double threshold;
	
	for (i = 0; i < NB_LOOPS; i++) {
		threshold = 0.01 * genrand();
		A = random_vector_dont_mangle_2(1.0);
		
		gradeBitmap1 = gradeBitmap(A, threshold);
		
		gmvA = new mv(A);
		
		gradeBitmap2 = gradeBitmap(gmvA, threshold);
		
		// check if grade bitmaps match
		if (gradeBitmap1 != gradeBitmap2) {
			Console.WriteLine("gradeBitmap() test failed (grade bitmap " + gradeBitmap1 + " vs " + gradeBitmap2 + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_gradeBitmap_dont_mangle_484(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	mv gmvA;
	bivector A;
	int i;
	int gradeBitmap1, gradeBitmap2;
	double threshold;
	
	for (i = 0; i < NB_LOOPS; i++) {
		threshold = 0.01 * genrand();
		A = random_bivector_dont_mangle_4(1.0);
		
		gradeBitmap1 = gradeBitmap(A, threshold);
		
		gmvA = new mv(A);
		
		gradeBitmap2 = gradeBitmap(gmvA, threshold);
		
		// check if grade bitmaps match
		if (gradeBitmap1 != gradeBitmap2) {
			Console.WriteLine("gradeBitmap() test failed (grade bitmap " + gradeBitmap1 + " vs " + gradeBitmap2 + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_gradeBitmap_dont_mangle_485(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	mv gmvA;
	trivector A;
	int i;
	int gradeBitmap1, gradeBitmap2;
	double threshold;
	
	for (i = 0; i < NB_LOOPS; i++) {
		threshold = 0.01 * genrand();
		A = random_trivector_dont_mangle_7(1.0);
		
		gradeBitmap1 = gradeBitmap(A, threshold);
		
		gmvA = new mv(A);
		
		gradeBitmap2 = gradeBitmap(gmvA, threshold);
		
		// check if grade bitmaps match
		if (gradeBitmap1 != gradeBitmap2) {
			Console.WriteLine("gradeBitmap() test failed (grade bitmap " + gradeBitmap1 + " vs " + gradeBitmap2 + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_gradeBitmap_dont_mangle_486(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	mv gmvA;
	e1_t A;
	int i;
	int gradeBitmap1, gradeBitmap2;
	double threshold;
	
	for (i = 0; i < NB_LOOPS; i++) {
		threshold = 0.01 * genrand();
		A = random_e1_t_dont_mangle_10(1.0);
		
		gradeBitmap1 = gradeBitmap(A, threshold);
		
		gmvA = new mv(A);
		
		gradeBitmap2 = gradeBitmap(gmvA, threshold);
		
		// check if grade bitmaps match
		if (gradeBitmap1 != gradeBitmap2) {
			Console.WriteLine("gradeBitmap() test failed (grade bitmap " + gradeBitmap1 + " vs " + gradeBitmap2 + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_gradeBitmap_dont_mangle_487(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	mv gmvA;
	e2_t A;
	int i;
	int gradeBitmap1, gradeBitmap2;
	double threshold;
	
	for (i = 0; i < NB_LOOPS; i++) {
		threshold = 0.01 * genrand();
		A = random_e2_t_dont_mangle_11(1.0);
		
		gradeBitmap1 = gradeBitmap(A, threshold);
		
		gmvA = new mv(A);
		
		gradeBitmap2 = gradeBitmap(gmvA, threshold);
		
		// check if grade bitmaps match
		if (gradeBitmap1 != gradeBitmap2) {
			Console.WriteLine("gradeBitmap() test failed (grade bitmap " + gradeBitmap1 + " vs " + gradeBitmap2 + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_gradeBitmap_dont_mangle_488(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	mv gmvA;
	I3_t A;
	int i;
	int gradeBitmap1, gradeBitmap2;
	double threshold;
	
	for (i = 0; i < NB_LOOPS; i++) {
		threshold = 0.01 * genrand();
		A = random_I3_t_dont_mangle_13(1.0);
		
		gradeBitmap1 = gradeBitmap(A, threshold);
		
		gmvA = new mv(A);
		
		gradeBitmap2 = gradeBitmap(gmvA, threshold);
		
		// check if grade bitmaps match
		if (gradeBitmap1 != gradeBitmap2) {
			Console.WriteLine("gradeBitmap() test failed (grade bitmap " + gradeBitmap1 + " vs " + gradeBitmap2 + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_hp_dont_mangle_489(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, C, D, dif;
	int i, g;
	int basisVectorBitmap = -1;
	double s;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		s = genrand();
		g = (int)(genrand() * 3.5);
		A = random_versor_dont_mangle_1_returns_mv(s, g, basisVectorBitmap);
		
		// copy it to another versor
		B = new mv(A);
		
		// set coordinates of B to random values (which may not be zero)
		for (g = 0; g < 4; g++) {
			if (B.m_c[g] == null) continue;
			for (int e = 0; e < B.m_c[g].Length; e++) {
				B.m_c[g][e] = 0.5 + genrand();
			}
		}
		
		// do hadamard product
		C = hp(A, B);
		
		// invert coordinates of B manually
		for (g = 0; g < 4; g++) {
			if (B.m_c[g] == null) continue;
			for (int e = 0; e < B.m_c[g].Length; e++) {
				B.m_c[g][e] = 1.0 / B.m_c[g][e];
			}
		}

		// do inverse hadamard product
		D = hp(C, B);
		
		// check if equal to original:
		dif = subtract(A, D);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("hp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_hp_dont_mangle_490(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	mv gmvA, gmvB, gmvC, gmvD;
	vector A;
	vector B;
	vector C;
	int i;
	mv dif;


	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_vector_dont_mangle_2(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		
		// do hadamard product (SMV)
		C = hp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = hp(gmvA, gmvB);
		
		// check if equal to original:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("hp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_hp_dont_mangle_491(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	mv gmvA, gmvB, gmvC, gmvD;
	bivector A;
	bivector B;
	bivector C;
	int i;
	double vC, vD, Q;


	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_bivector_dont_mangle_4(1.0);
		B = random_bivector_dont_mangle_4(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		
		// do hadamard product (SMV)
		C = hp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = hp(gmvA, gmvB);
		
		vC = gmvC.get_scalar();
		vD = gmvD.get_scalar();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1();
		vD = gmvD.get_e1();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e2();
		vD = gmvD.get_e2();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e3();
		vD = gmvD.get_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e2();
		vD = gmvD.get_e1_e2();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e3();
		vD = gmvD.get_e1_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e2_e3();
		vD = gmvD.get_e2_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e2_e3();
		vD = gmvD.get_e1_e2_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
	}
	return 1; // success
}

static int test_hp_dont_mangle_492(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	mv gmvA, gmvB, gmvC, gmvD;
	rotor A;
	rotor B;
	rotor C;
	int i;
	double vC, vD, Q;


	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_rotor_dont_mangle_8(1.0);
		B = random_rotor_dont_mangle_8(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		
		// do hadamard product (SMV)
		C = hp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = hp(gmvA, gmvB);
		
		vC = gmvC.get_scalar();
		vD = gmvD.get_scalar();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1();
		vD = gmvD.get_e1();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e2();
		vD = gmvD.get_e2();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e3();
		vD = gmvD.get_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e2();
		vD = gmvD.get_e1_e2();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e3();
		vD = gmvD.get_e1_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e2_e3();
		vD = gmvD.get_e2_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e2_e3();
		vD = gmvD.get_e1_e2_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
	}
	return 1; // success
}

static int test_hp_dont_mangle_493(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	mv gmvA, gmvB, gmvC, gmvD;
	bivector A;
	rotor B;
	bivector C;
	int i;
	double vC, vD, Q;


	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_bivector_dont_mangle_4(1.0);
		B = random_rotor_dont_mangle_8(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		
		// do hadamard product (SMV)
		C = hp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = hp(gmvA, gmvB);
		
		vC = gmvC.get_scalar();
		vD = gmvD.get_scalar();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1();
		vD = gmvD.get_e1();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e2();
		vD = gmvD.get_e2();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e3();
		vD = gmvD.get_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e2();
		vD = gmvD.get_e1_e2();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e3();
		vD = gmvD.get_e1_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e2_e3();
		vD = gmvD.get_e2_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e2_e3();
		vD = gmvD.get_e1_e2_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
	}
	return 1; // success
}

static int test_hp_dont_mangle_494(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	mv gmvA, gmvB, gmvC, gmvD;
	trivector A;
	trivector B;
	trivector C;
	int i;
	mv dif;


	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_trivector_dont_mangle_7(1.0);
		B = random_trivector_dont_mangle_7(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		
		// do hadamard product (SMV)
		C = hp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = hp(gmvA, gmvB);
		
		// check if equal to original:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("hp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_hp_dont_mangle_495(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	mv gmvA, gmvB, gmvC, gmvD;
	trivector A;
	oddVersor B;
	trivector C;
	int i;
	mv dif;


	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_trivector_dont_mangle_7(1.0);
		B = random_oddVersor_dont_mangle_93(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		
		// do hadamard product (SMV)
		C = hp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = hp(gmvA, gmvB);
		
		// check if equal to original:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("hp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_hp_dont_mangle_496(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	mv gmvA, gmvB, gmvC, gmvD;
	rotor A;
	bivector B;
	bivector C;
	int i;
	double vC, vD, Q;


	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_rotor_dont_mangle_8(1.0);
		B = random_bivector_dont_mangle_4(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		
		// do hadamard product (SMV)
		C = hp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = hp(gmvA, gmvB);
		
		vC = gmvC.get_scalar();
		vD = gmvD.get_scalar();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1();
		vD = gmvD.get_e1();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e2();
		vD = gmvD.get_e2();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e3();
		vD = gmvD.get_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e2();
		vD = gmvD.get_e1_e2();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e3();
		vD = gmvD.get_e1_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e2_e3();
		vD = gmvD.get_e2_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e2_e3();
		vD = gmvD.get_e1_e2_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("hp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
	}
	return 1; // success
}

static int test_hp_dont_mangle_497(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	mv gmvA, gmvB, gmvC, gmvD;
	e1_t A;
	e1_t B;
	e1_t C;
	int i;
	mv dif;


	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_e1_t_dont_mangle_10(1.0);
		B = random_e1_t_dont_mangle_10(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		
		// do hadamard product (SMV)
		C = hp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = hp(gmvA, gmvB);
		
		// check if equal to original:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("hp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_hp_dont_mangle_498(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	mv gmvA, gmvB, gmvC, gmvD;
	e2_t A;
	e3_t B;
	double C;
	int i;
	mv dif;


	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_e2_t_dont_mangle_11(1.0);
		B = random_e3_t_dont_mangle_75(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		
		// do hadamard product (SMV)
		C = hp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = hp(gmvA, gmvB);
		
		// check if equal to original:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("hp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_hp_dont_mangle_499(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	mv gmvA, gmvB, gmvC, gmvD;
	oddVersor A;
	I3_t B;
	trivector C;
	int i;
	mv dif;


	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_oddVersor_dont_mangle_93(1.0);
		B = random_I3_t_dont_mangle_13(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		
		// do hadamard product (SMV)
		C = hp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = hp(gmvA, gmvB);
		
		// check if equal to original:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("hp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_ihp_dont_mangle_500(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, C, D, dif;
	int i, g;
	int basisVectorBitmap = -1;
	double s;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		s = genrand();
		g = (int)(genrand() * 3.5);
		A = random_versor_dont_mangle_1_returns_mv(s, g, basisVectorBitmap);
		
		// copy it to another versor
		B = new mv(A);
		
		// set coordinates of B to random values (which may not be zero)
		for (g = 0; g < 4; g++) {
			if (B.m_c[g] == null) continue;
			for (int e = 0; e < B.m_c[g].Length; e++) {
				B.m_c[g][e] = 0.5 + genrand();
			}
		}
		
		// do hadamard product
		C = ihp(A, B);
		
		// invert coordinates of B manually
		for (g = 0; g < 4; g++) {
			if (B.m_c[g] == null) continue;
			for (int e = 0; e < B.m_c[g].Length; e++) {
				B.m_c[g][e] = 1.0 / B.m_c[g][e];
			}
		}

		// do inverse hadamard product
		D = ihp(C, B);
		
		// check if equal to original:
		dif = subtract(A, D);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("ihp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_ihp_dont_mangle_501(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	mv gmvA, gmvB, gmvC, gmvD;
	vector A;
	vector B;
	vector C;
	int i;
	bool ok;
	mv dif;

	// initialize a lookup table which tells us which coordinates should not be non-null	
	bool[] nonNull2 = new bool[8];
	for (i = 0; i < 100; i++) {
		B = random_vector_dont_mangle_2(1.0);
		gmvB = new mv(B);
		int nonNull2idx = 0;
		for (int g = 0; g < 4; g++) {
			if (gmvB.m_c[g] == null) continue;
			for (int e = 0; e < gmvB.m_c[g].Length; e++) {
				if (gmvB.m_c[g][e] != 0.0) nonNull2[nonNull2idx] = true;
				nonNull2idx++;
			}
		}
	}

	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_vector_dont_mangle_2(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		// make sure that none of the 'B' coordinates are 0, because then inverse hadamard won't work
		ok = true;

		{		
			int nonNull2idx = 0;
			for (int g = 0; g < 4; g++) {
				if (gmvB.m_c[g] == null) continue;
				for (int e = 0; e < gmvB.m_c[g].Length; e++) {
					if (gmvB.m_c[g][e] == 0.0) {
						if (nonNull2[nonNull2idx]) ok = false;
						else gmvB.m_c[g][e] = 1E+306;
						nonNull2idx++;
					}
				}
			}
		}
		
		if (!ok) continue; // some SMV coordinate was set to 0 by the random generator (very rare)
		
		// do hadamard product (SMV)
		C = ihp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = ihp(gmvA, gmvB);
		
		// check if equal to original:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("ihp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_ihp_dont_mangle_502(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	mv gmvA, gmvB, gmvC, gmvD;
	bivector A;
	bivector B;
	bivector C;
	int i;
	bool ok;
	double vC, vD, Q;

	// initialize a lookup table which tells us which coordinates should not be non-null	
	bool[] nonNull2 = new bool[8];
	for (i = 0; i < 100; i++) {
		B = random_bivector_dont_mangle_4(1.0);
		gmvB = new mv(B);
		int nonNull2idx = 0;
		for (int g = 0; g < 4; g++) {
			if (gmvB.m_c[g] == null) continue;
			for (int e = 0; e < gmvB.m_c[g].Length; e++) {
				if (gmvB.m_c[g][e] != 0.0) nonNull2[nonNull2idx] = true;
				nonNull2idx++;
			}
		}
	}

	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_bivector_dont_mangle_4(1.0);
		B = random_bivector_dont_mangle_4(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		// make sure that none of the 'B' coordinates are 0, because then inverse hadamard won't work
		ok = true;

		{		
			int nonNull2idx = 0;
			for (int g = 0; g < 4; g++) {
				if (gmvB.m_c[g] == null) continue;
				for (int e = 0; e < gmvB.m_c[g].Length; e++) {
					if (gmvB.m_c[g][e] == 0.0) {
						if (nonNull2[nonNull2idx]) ok = false;
						else gmvB.m_c[g][e] = 1E+306;
						nonNull2idx++;
					}
				}
			}
		}
		
		if (!ok) continue; // some SMV coordinate was set to 0 by the random generator (very rare)
		
		// do hadamard product (SMV)
		C = ihp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = ihp(gmvA, gmvB);
		
		vC = gmvC.get_scalar();
		vD = gmvD.get_scalar();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1();
		vD = gmvD.get_e1();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e2();
		vD = gmvD.get_e2();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e3();
		vD = gmvD.get_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e2();
		vD = gmvD.get_e1_e2();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e3();
		vD = gmvD.get_e1_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e2_e3();
		vD = gmvD.get_e2_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e2_e3();
		vD = gmvD.get_e1_e2_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
	}
	return 1; // success
}

static int test_ihp_dont_mangle_503(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	mv gmvA, gmvB, gmvC, gmvD;
	rotor A;
	rotor B;
	rotor C;
	int i;
	bool ok;
	double vC, vD, Q;

	// initialize a lookup table which tells us which coordinates should not be non-null	
	bool[] nonNull2 = new bool[8];
	for (i = 0; i < 100; i++) {
		B = random_rotor_dont_mangle_8(1.0);
		gmvB = new mv(B);
		int nonNull2idx = 0;
		for (int g = 0; g < 4; g++) {
			if (gmvB.m_c[g] == null) continue;
			for (int e = 0; e < gmvB.m_c[g].Length; e++) {
				if (gmvB.m_c[g][e] != 0.0) nonNull2[nonNull2idx] = true;
				nonNull2idx++;
			}
		}
	}

	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_rotor_dont_mangle_8(1.0);
		B = random_rotor_dont_mangle_8(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		// make sure that none of the 'B' coordinates are 0, because then inverse hadamard won't work
		ok = true;

		{		
			int nonNull2idx = 0;
			for (int g = 0; g < 4; g++) {
				if (gmvB.m_c[g] == null) continue;
				for (int e = 0; e < gmvB.m_c[g].Length; e++) {
					if (gmvB.m_c[g][e] == 0.0) {
						if (nonNull2[nonNull2idx]) ok = false;
						else gmvB.m_c[g][e] = 1E+306;
						nonNull2idx++;
					}
				}
			}
		}
		
		if (!ok) continue; // some SMV coordinate was set to 0 by the random generator (very rare)
		
		// do hadamard product (SMV)
		C = ihp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = ihp(gmvA, gmvB);
		
		vC = gmvC.get_scalar();
		vD = gmvD.get_scalar();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1();
		vD = gmvD.get_e1();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e2();
		vD = gmvD.get_e2();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e3();
		vD = gmvD.get_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e2();
		vD = gmvD.get_e1_e2();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e3();
		vD = gmvD.get_e1_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e2_e3();
		vD = gmvD.get_e2_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e2_e3();
		vD = gmvD.get_e1_e2_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
	}
	return 1; // success
}

static int test_ihp_dont_mangle_504(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	mv gmvA, gmvB, gmvC, gmvD;
	bivector A;
	rotor B;
	bivector C;
	int i;
	bool ok;
	double vC, vD, Q;

	// initialize a lookup table which tells us which coordinates should not be non-null	
	bool[] nonNull2 = new bool[8];
	for (i = 0; i < 100; i++) {
		B = random_rotor_dont_mangle_8(1.0);
		gmvB = new mv(B);
		int nonNull2idx = 0;
		for (int g = 0; g < 4; g++) {
			if (gmvB.m_c[g] == null) continue;
			for (int e = 0; e < gmvB.m_c[g].Length; e++) {
				if (gmvB.m_c[g][e] != 0.0) nonNull2[nonNull2idx] = true;
				nonNull2idx++;
			}
		}
	}

	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_bivector_dont_mangle_4(1.0);
		B = random_rotor_dont_mangle_8(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		// make sure that none of the 'B' coordinates are 0, because then inverse hadamard won't work
		ok = true;

		{		
			int nonNull2idx = 0;
			for (int g = 0; g < 4; g++) {
				if (gmvB.m_c[g] == null) continue;
				for (int e = 0; e < gmvB.m_c[g].Length; e++) {
					if (gmvB.m_c[g][e] == 0.0) {
						if (nonNull2[nonNull2idx]) ok = false;
						else gmvB.m_c[g][e] = 1E+306;
						nonNull2idx++;
					}
				}
			}
		}
		
		if (!ok) continue; // some SMV coordinate was set to 0 by the random generator (very rare)
		
		// do hadamard product (SMV)
		C = ihp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = ihp(gmvA, gmvB);
		
		vC = gmvC.get_scalar();
		vD = gmvD.get_scalar();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1();
		vD = gmvD.get_e1();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e2();
		vD = gmvD.get_e2();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e3();
		vD = gmvD.get_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e2();
		vD = gmvD.get_e1_e2();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e3();
		vD = gmvD.get_e1_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e2_e3();
		vD = gmvD.get_e2_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e2_e3();
		vD = gmvD.get_e1_e2_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
	}
	return 1; // success
}

static int test_ihp_dont_mangle_505(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	mv gmvA, gmvB, gmvC, gmvD;
	rotor A;
	bivector B;
	bivector C;
	int i;
	bool ok;
	double vC, vD, Q;

	// initialize a lookup table which tells us which coordinates should not be non-null	
	bool[] nonNull2 = new bool[8];
	for (i = 0; i < 100; i++) {
		B = random_bivector_dont_mangle_4(1.0);
		gmvB = new mv(B);
		int nonNull2idx = 0;
		for (int g = 0; g < 4; g++) {
			if (gmvB.m_c[g] == null) continue;
			for (int e = 0; e < gmvB.m_c[g].Length; e++) {
				if (gmvB.m_c[g][e] != 0.0) nonNull2[nonNull2idx] = true;
				nonNull2idx++;
			}
		}
	}

	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_rotor_dont_mangle_8(1.0);
		B = random_bivector_dont_mangle_4(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		// make sure that none of the 'B' coordinates are 0, because then inverse hadamard won't work
		ok = true;

		{		
			int nonNull2idx = 0;
			for (int g = 0; g < 4; g++) {
				if (gmvB.m_c[g] == null) continue;
				for (int e = 0; e < gmvB.m_c[g].Length; e++) {
					if (gmvB.m_c[g][e] == 0.0) {
						if (nonNull2[nonNull2idx]) ok = false;
						else gmvB.m_c[g][e] = 1E+306;
						nonNull2idx++;
					}
				}
			}
		}
		
		if (!ok) continue; // some SMV coordinate was set to 0 by the random generator (very rare)
		
		// do hadamard product (SMV)
		C = ihp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = ihp(gmvA, gmvB);
		
		vC = gmvC.get_scalar();
		vD = gmvD.get_scalar();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1();
		vD = gmvD.get_e1();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e2();
		vD = gmvD.get_e2();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e3();
		vD = gmvD.get_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e2();
		vD = gmvD.get_e1_e2();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e3();
		vD = gmvD.get_e1_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e2_e3();
		vD = gmvD.get_e2_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
		vC = gmvC.get_e1_e2_e3();
		vD = gmvD.get_e1_e2_e3();
		if ((vD < -1E-14) || (vD > 1E-14)) {
			Q = vC / vD;
			if (Q < 0.0) Q = -Q;
			if (((Q - 1.0) < -1E-14) || ((Q - 1.0) > 1E-14)) {
				Console.WriteLine("ihp() test failed (using alternate test, Q = " + (double)Q + ")\n");
				return 0; // failure
			}
		}
	}
	return 1; // success
}

static int test_ihp_dont_mangle_506(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	mv gmvA, gmvB, gmvC, gmvD;
	trivector A;
	oddVersor B;
	trivector C;
	int i;
	bool ok;
	mv dif;

	// initialize a lookup table which tells us which coordinates should not be non-null	
	bool[] nonNull2 = new bool[8];
	for (i = 0; i < 100; i++) {
		B = random_oddVersor_dont_mangle_93(1.0);
		gmvB = new mv(B);
		int nonNull2idx = 0;
		for (int g = 0; g < 4; g++) {
			if (gmvB.m_c[g] == null) continue;
			for (int e = 0; e < gmvB.m_c[g].Length; e++) {
				if (gmvB.m_c[g][e] != 0.0) nonNull2[nonNull2idx] = true;
				nonNull2idx++;
			}
		}
	}

	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_trivector_dont_mangle_7(1.0);
		B = random_oddVersor_dont_mangle_93(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		// make sure that none of the 'B' coordinates are 0, because then inverse hadamard won't work
		ok = true;

		{		
			int nonNull2idx = 0;
			for (int g = 0; g < 4; g++) {
				if (gmvB.m_c[g] == null) continue;
				for (int e = 0; e < gmvB.m_c[g].Length; e++) {
					if (gmvB.m_c[g][e] == 0.0) {
						if (nonNull2[nonNull2idx]) ok = false;
						else gmvB.m_c[g][e] = 1E+306;
						nonNull2idx++;
					}
				}
			}
		}
		
		if (!ok) continue; // some SMV coordinate was set to 0 by the random generator (very rare)
		
		// do hadamard product (SMV)
		C = ihp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = ihp(gmvA, gmvB);
		
		// check if equal to original:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("ihp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_ihp_dont_mangle_507(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	mv gmvA, gmvB, gmvC, gmvD;
	vector A;
	e1_t B;
	vector C;
	int i;
	bool ok;
	mv dif;

	// initialize a lookup table which tells us which coordinates should not be non-null	
	bool[] nonNull2 = new bool[8];
	for (i = 0; i < 100; i++) {
		B = random_e1_t_dont_mangle_10(1.0);
		gmvB = new mv(B);
		int nonNull2idx = 0;
		for (int g = 0; g < 4; g++) {
			if (gmvB.m_c[g] == null) continue;
			for (int e = 0; e < gmvB.m_c[g].Length; e++) {
				if (gmvB.m_c[g][e] != 0.0) nonNull2[nonNull2idx] = true;
				nonNull2idx++;
			}
		}
	}

	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_e1_t_dont_mangle_10(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		// make sure that none of the 'B' coordinates are 0, because then inverse hadamard won't work
		ok = true;

		{		
			int nonNull2idx = 0;
			for (int g = 0; g < 4; g++) {
				if (gmvB.m_c[g] == null) continue;
				for (int e = 0; e < gmvB.m_c[g].Length; e++) {
					if (gmvB.m_c[g][e] == 0.0) {
						if (nonNull2[nonNull2idx]) ok = false;
						else gmvB.m_c[g][e] = 1E+306;
						nonNull2idx++;
					}
				}
			}
		}
		
		if (!ok) continue; // some SMV coordinate was set to 0 by the random generator (very rare)
		
		// do hadamard product (SMV)
		C = ihp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = ihp(gmvA, gmvB);
		
		// check if equal to original:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("ihp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_ihp_dont_mangle_508(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	mv gmvA, gmvB, gmvC, gmvD;
	e2_t A;
	e3_t B;
	double C;
	int i;
	bool ok;
	mv dif;

	// initialize a lookup table which tells us which coordinates should not be non-null	
	bool[] nonNull2 = new bool[8];
	for (i = 0; i < 100; i++) {
		B = random_e3_t_dont_mangle_75(1.0);
		gmvB = new mv(B);
		int nonNull2idx = 0;
		for (int g = 0; g < 4; g++) {
			if (gmvB.m_c[g] == null) continue;
			for (int e = 0; e < gmvB.m_c[g].Length; e++) {
				if (gmvB.m_c[g][e] != 0.0) nonNull2[nonNull2idx] = true;
				nonNull2idx++;
			}
		}
	}

	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_e2_t_dont_mangle_11(1.0);
		B = random_e3_t_dont_mangle_75(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		// make sure that none of the 'B' coordinates are 0, because then inverse hadamard won't work
		ok = true;

		{		
			int nonNull2idx = 0;
			for (int g = 0; g < 4; g++) {
				if (gmvB.m_c[g] == null) continue;
				for (int e = 0; e < gmvB.m_c[g].Length; e++) {
					if (gmvB.m_c[g][e] == 0.0) {
						if (nonNull2[nonNull2idx]) ok = false;
						else gmvB.m_c[g][e] = 1E+306;
						nonNull2idx++;
					}
				}
			}
		}
		
		if (!ok) continue; // some SMV coordinate was set to 0 by the random generator (very rare)
		
		// do hadamard product (SMV)
		C = ihp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = ihp(gmvA, gmvB);
		
		// check if equal to original:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("ihp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_ihp_dont_mangle_509(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	mv gmvA, gmvB, gmvC, gmvD;
	trivector A;
	I3_t B;
	trivector C;
	int i;
	bool ok;
	mv dif;

	// initialize a lookup table which tells us which coordinates should not be non-null	
	bool[] nonNull2 = new bool[8];
	for (i = 0; i < 100; i++) {
		B = random_I3_t_dont_mangle_13(1.0);
		gmvB = new mv(B);
		int nonNull2idx = 0;
		for (int g = 0; g < 4; g++) {
			if (gmvB.m_c[g] == null) continue;
			for (int e = 0; e < gmvB.m_c[g].Length; e++) {
				if (gmvB.m_c[g][e] != 0.0) nonNull2[nonNull2idx] = true;
				nonNull2idx++;
			}
		}
	}

	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_trivector_dont_mangle_7(1.0);
		B = random_I3_t_dont_mangle_13(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
		
		// make sure that none of the 'B' coordinates are 0, because then inverse hadamard won't work
		ok = true;

		{		
			int nonNull2idx = 0;
			for (int g = 0; g < 4; g++) {
				if (gmvB.m_c[g] == null) continue;
				for (int e = 0; e < gmvB.m_c[g].Length; e++) {
					if (gmvB.m_c[g][e] == 0.0) {
						if (nonNull2[nonNull2idx]) ok = false;
						else gmvB.m_c[g][e] = 1E+306;
						nonNull2idx++;
					}
				}
			}
		}
		
		if (!ok) continue; // some SMV coordinate was set to 0 by the random generator (very rare)
		
		// do hadamard product (SMV)
		C = ihp(A, B);
		gmvC = new mv(C);
		
		// do hadamard product (GMV)
		gmvD = ihp(gmvA, gmvB);
		
		// check if equal to original:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("ihp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_increment_dont_mangle_510(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, C, D, one;
	int i, g;
	int basisVectorBitmap = -1;

	one = new mv(1.0);

	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		g = (int)(genrand() * 3.5);
		A = random_versor_dont_mangle_1_returns_mv(genrand() + 0.5, g, basisVectorBitmap);
		
		// increment/decrement
		B = increment(A);
		
		// undo the increment/decrement
		C = subtract(B, one);
		
		// see if (A - (B-1)) == 0
		D = subtract(A, C);
		
		if (D.LargestCoordinate() > 1E-14) {
			Console.WriteLine("increment() test failed (largestCoordinate of D = " + (double)D.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_increment_dont_mangle_511(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	rotor C;
	mv gA, gC1, gC2;
	double s;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_bivector_dont_mangle_4(s);
		
		// increment or decrement
		C = increment(A);
		
		// convert A and C to gmv
		gA = new mv(A);
		gC1 = new mv(C);
		
		gC2 = increment(gA);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-14) {
			Console.WriteLine("increment() test failed\n");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_increment_dont_mangle_512(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	rotor C;
	mv gA, gC1, gC2;
	double s;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_rotor_dont_mangle_8(s);
		
		// increment or decrement
		C = increment(A);
		
		// convert A and C to gmv
		gA = new mv(A);
		gC1 = new mv(C);
		
		gC2 = increment(gA);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-14) {
			Console.WriteLine("increment() test failed\n");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_decrement_dont_mangle_513(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, C, D, one;
	int i, g;
	int basisVectorBitmap = -1;

	one = new mv(-1.0);

	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		g = (int)(genrand() * 3.5);
		A = random_versor_dont_mangle_1_returns_mv(genrand() + 0.5, g, basisVectorBitmap);
		
		// increment/decrement
		B = decrement(A);
		
		// undo the increment/decrement
		C = subtract(B, one);
		
		// see if (A - (B-1)) == 0
		D = subtract(A, C);
		
		if (D.LargestCoordinate() > 1E-14) {
			Console.WriteLine("decrement() test failed (largestCoordinate of D = " + (double)D.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_decrement_dont_mangle_514(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	rotor C;
	mv gA, gC1, gC2;
	double s;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_bivector_dont_mangle_4(s);
		
		// increment or decrement
		C = decrement(A);
		
		// convert A and C to gmv
		gA = new mv(A);
		gC1 = new mv(C);
		
		gC2 = decrement(gA);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-14) {
			Console.WriteLine("decrement() test failed\n");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_decrement_dont_mangle_515(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	rotor C;
	mv gA, gC1, gC2;
	double s;
	
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_rotor_dont_mangle_8(s);
		
		// increment or decrement
		C = decrement(A);
		
		// convert A and C to gmv
		gA = new mv(A);
		gC1 = new mv(C);
		
		gC2 = decrement(gA);
		
		// see if result is equal up to precision:
		gA = subtract(gC1, gC2);
		if (gA.LargestCoordinate() > 1E-14) {
			Console.WriteLine("decrement() test failed\n");
			return 0; // failure
		}		
	}
	return 1; // success
}

static int test_sp_dont_mangle_516(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, C, D, E, dif;
	int i, ga, gb, gd;
	double s;
	int basisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		ga = (int)(genrand() * 3.5);
		A = random_blade_dont_mangle_23_returns_mv(s, ga, basisVectorBitmap);
		
		s = genrand();
		gb = (int)(genrand() * 3.5);
		B = random_blade_dont_mangle_23_returns_mv(s, gb, basisVectorBitmap);
		
		// compute product using sp()
		C = new mv(sp(A, B));
		
		// simulate product using geometric product & grade part selection
		D = gp(A, B);
		gd = (ga > gb) ? ga - gb : gb - ga;
		E = extractGrade(D, Grades[0]);

		// check if equal:
		dif = subtract(C, E);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("sp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_lc_dont_mangle_517(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, C, D, E, dif;
	int i, ga, gb, gd;
	double s;
	int basisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		ga = (int)(genrand() * 3.5);
		A = random_blade_dont_mangle_23_returns_mv(s, ga, basisVectorBitmap);
		
		s = genrand();
		gb = (int)(genrand() * 3.5);
		B = random_blade_dont_mangle_23_returns_mv(s, gb, basisVectorBitmap);
		
		// compute product using lc()
		C = new mv(lc(A, B));
		
		// simulate product using geometric product & grade part selection
		D = gp(A, B);
		gd = (ga > gb) ? ga - gb : gb - ga;
		if (ga > gb) E = new mv(0.0);
		else E = extractGrade(D, Grades[gd]);

		// check if equal:
		dif = subtract(C, E);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("lc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_rc_dont_mangle_518(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, C, D, E, dif;
	int i, ga, gb, gd;
	double s;
	int basisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		ga = (int)(genrand() * 3.5);
		A = random_blade_dont_mangle_23_returns_mv(s, ga, basisVectorBitmap);
		
		s = genrand();
		gb = (int)(genrand() * 3.5);
		B = random_blade_dont_mangle_23_returns_mv(s, gb, basisVectorBitmap);
		
		// compute product using rc()
		C = new mv(rc(A, B));
		
		// simulate product using geometric product & grade part selection
		D = gp(A, B);
		gd = (ga > gb) ? ga - gb : gb - ga;
		if (ga < gb) E = new mv(0.0);
		else E = extractGrade(D, Grades[gd]);

		// check if equal:
		dif = subtract(C, E);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("rc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_hip_dont_mangle_519(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, C, D, E, dif;
	int i, ga, gb, gd;
	double s;
	int basisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		ga = (int)(genrand() * 3.5);
		A = random_blade_dont_mangle_23_returns_mv(s, ga, basisVectorBitmap);
		
		s = genrand();
		gb = (int)(genrand() * 3.5);
		B = random_blade_dont_mangle_23_returns_mv(s, gb, basisVectorBitmap);
		
		// compute product using hip()
		C = new mv(hip(A, B));
		
		// simulate product using geometric product & grade part selection
		D = gp(A, B);
		gd = (ga > gb) ? ga - gb : gb - ga;
		if ((ga == 0) || (gb == 0)) E = new mv(0.0);
		else E = extractGrade(D, Grades[gd]);

		// check if equal:
		dif = subtract(C, E);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("hip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_mhip_dont_mangle_520(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, C, D, E, dif;
	int i, ga, gb, gd;
	double s;
	int basisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		ga = (int)(genrand() * 3.5);
		A = random_blade_dont_mangle_23_returns_mv(s, ga, basisVectorBitmap);
		
		s = genrand();
		gb = (int)(genrand() * 3.5);
		B = random_blade_dont_mangle_23_returns_mv(s, gb, basisVectorBitmap);
		
		// compute product using mhip()
		C = new mv(mhip(A, B));
		
		// simulate product using geometric product & grade part selection
		D = gp(A, B);
		gd = (ga > gb) ? ga - gb : gb - ga;
		E = extractGrade(D, Grades[gd]);

		// check if equal:
		dif = subtract(C, E);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("mhip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_sp_dont_mangle_521(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	vector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_vector_dont_mangle_2(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using sp()
		gmvC = new mv(sp(A, B));
		
		// perform GMV version 
		gmvD = new mv(sp(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("sp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_lc_dont_mangle_522(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	vector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_vector_dont_mangle_2(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using lc()
		gmvC = new mv(lc(A, B));
		
		// perform GMV version 
		gmvD = new mv(lc(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("lc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_rc_dont_mangle_523(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	vector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_vector_dont_mangle_2(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using rc()
		gmvC = new mv(rc(A, B));
		
		// perform GMV version 
		gmvD = new mv(rc(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("rc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_hip_dont_mangle_524(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	vector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_vector_dont_mangle_2(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using hip()
		gmvC = new mv(hip(A, B));
		
		// perform GMV version 
		gmvD = new mv(hip(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("hip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_mhip_dont_mangle_525(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	vector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_vector_dont_mangle_2(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using mhip()
		gmvC = new mv(mhip(A, B));
		
		// perform GMV version 
		gmvD = new mv(mhip(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("mhip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_sp_dont_mangle_526(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	vector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_bivector_dont_mangle_4(1.0);
		B = random_vector_dont_mangle_2(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using sp()
		gmvC = new mv(sp(A, B));
		
		// perform GMV version 
		gmvD = new mv(sp(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("sp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_lc_dont_mangle_527(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	vector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_bivector_dont_mangle_4(1.0);
		B = random_vector_dont_mangle_2(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using lc()
		gmvC = new mv(lc(A, B));
		
		// perform GMV version 
		gmvD = new mv(lc(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("lc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_rc_dont_mangle_528(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	vector B;
//	vector C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_bivector_dont_mangle_4(1.0);
		B = random_vector_dont_mangle_2(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using rc()
		gmvC = new mv(rc(A, B));
		
		// perform GMV version 
		gmvD = new mv(rc(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("rc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_hip_dont_mangle_529(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	vector B;
//	vector C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_bivector_dont_mangle_4(1.0);
		B = random_vector_dont_mangle_2(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using hip()
		gmvC = new mv(hip(A, B));
		
		// perform GMV version 
		gmvD = new mv(hip(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("hip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_mhip_dont_mangle_530(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	vector B;
//	vector C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_bivector_dont_mangle_4(1.0);
		B = random_vector_dont_mangle_2(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using mhip()
		gmvC = new mv(mhip(A, B));
		
		// perform GMV version 
		gmvD = new mv(mhip(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("mhip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_sp_dont_mangle_531(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	trivector A;
	trivector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_trivector_dont_mangle_7(1.0);
		B = random_trivector_dont_mangle_7(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using sp()
		gmvC = new mv(sp(A, B));
		
		// perform GMV version 
		gmvD = new mv(sp(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("sp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_lc_dont_mangle_532(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	trivector A;
	trivector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_trivector_dont_mangle_7(1.0);
		B = random_trivector_dont_mangle_7(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using lc()
		gmvC = new mv(lc(A, B));
		
		// perform GMV version 
		gmvD = new mv(lc(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("lc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_rc_dont_mangle_533(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	trivector A;
	trivector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_trivector_dont_mangle_7(1.0);
		B = random_trivector_dont_mangle_7(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using rc()
		gmvC = new mv(rc(A, B));
		
		// perform GMV version 
		gmvD = new mv(rc(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("rc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_hip_dont_mangle_534(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	trivector A;
	trivector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_trivector_dont_mangle_7(1.0);
		B = random_trivector_dont_mangle_7(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using hip()
		gmvC = new mv(hip(A, B));
		
		// perform GMV version 
		gmvD = new mv(hip(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("hip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_mhip_dont_mangle_535(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	trivector A;
	trivector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_trivector_dont_mangle_7(1.0);
		B = random_trivector_dont_mangle_7(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using mhip()
		gmvC = new mv(mhip(A, B));
		
		// perform GMV version 
		gmvD = new mv(mhip(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("mhip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_sp_dont_mangle_536(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	bivector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_bivector_dont_mangle_4(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using sp()
		gmvC = new mv(sp(A, B));
		
		// perform GMV version 
		gmvD = new mv(sp(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("sp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_lc_dont_mangle_537(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	bivector B;
//	vector C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_bivector_dont_mangle_4(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using lc()
		gmvC = new mv(lc(A, B));
		
		// perform GMV version 
		gmvD = new mv(lc(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("lc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_rc_dont_mangle_538(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	bivector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_bivector_dont_mangle_4(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using rc()
		gmvC = new mv(rc(A, B));
		
		// perform GMV version 
		gmvD = new mv(rc(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("rc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_hip_dont_mangle_539(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	bivector B;
//	vector C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_bivector_dont_mangle_4(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using hip()
		gmvC = new mv(hip(A, B));
		
		// perform GMV version 
		gmvD = new mv(hip(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("hip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_mhip_dont_mangle_540(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	bivector B;
//	vector C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_bivector_dont_mangle_4(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using mhip()
		gmvC = new mv(mhip(A, B));
		
		// perform GMV version 
		gmvD = new mv(mhip(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("mhip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_sp_dont_mangle_541(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	rotor B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_rotor_dont_mangle_8(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using sp()
		gmvC = new mv(sp(A, B));
		
		// perform GMV version 
		gmvD = new mv(sp(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("sp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_lc_dont_mangle_542(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	rotor B;
//	vector C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_rotor_dont_mangle_8(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using lc()
		gmvC = new mv(lc(A, B));
		
		// perform GMV version 
		gmvD = new mv(lc(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("lc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_rc_dont_mangle_543(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	rotor B;
//	vector C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_rotor_dont_mangle_8(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using rc()
		gmvC = new mv(rc(A, B));
		
		// perform GMV version 
		gmvD = new mv(rc(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("rc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_hip_dont_mangle_544(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	rotor B;
//	vector C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_rotor_dont_mangle_8(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using hip()
		gmvC = new mv(hip(A, B));
		
		// perform GMV version 
		gmvD = new mv(hip(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("hip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_mhip_dont_mangle_545(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	rotor B;
//	vector C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_rotor_dont_mangle_8(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using mhip()
		gmvC = new mv(mhip(A, B));
		
		// perform GMV version 
		gmvD = new mv(mhip(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("mhip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_sp_dont_mangle_546(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	bivector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_rotor_dont_mangle_8(1.0);
		B = random_bivector_dont_mangle_4(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using sp()
		gmvC = new mv(sp(A, B));
		
		// perform GMV version 
		gmvD = new mv(sp(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("sp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_lc_dont_mangle_547(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	bivector B;
//	rotor C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_rotor_dont_mangle_8(1.0);
		B = random_bivector_dont_mangle_4(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using lc()
		gmvC = new mv(lc(A, B));
		
		// perform GMV version 
		gmvD = new mv(lc(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("lc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_rc_dont_mangle_548(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	bivector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_rotor_dont_mangle_8(1.0);
		B = random_bivector_dont_mangle_4(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using rc()
		gmvC = new mv(rc(A, B));
		
		// perform GMV version 
		gmvD = new mv(rc(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("rc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_hip_dont_mangle_549(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	bivector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_rotor_dont_mangle_8(1.0);
		B = random_bivector_dont_mangle_4(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using hip()
		gmvC = new mv(hip(A, B));
		
		// perform GMV version 
		gmvD = new mv(hip(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("hip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_mhip_dont_mangle_550(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	bivector B;
//	rotor C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_rotor_dont_mangle_8(1.0);
		B = random_bivector_dont_mangle_4(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using mhip()
		gmvC = new mv(mhip(A, B));
		
		// perform GMV version 
		gmvD = new mv(mhip(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("mhip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_sp_dont_mangle_551(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	trivector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_rotor_dont_mangle_8(1.0);
		B = random_trivector_dont_mangle_7(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using sp()
		gmvC = new mv(sp(A, B));
		
		// perform GMV version 
		gmvD = new mv(sp(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("sp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_lc_dont_mangle_552(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	trivector B;
//	oddVersor C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_rotor_dont_mangle_8(1.0);
		B = random_trivector_dont_mangle_7(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using lc()
		gmvC = new mv(lc(A, B));
		
		// perform GMV version 
		gmvD = new mv(lc(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("lc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_rc_dont_mangle_553(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	trivector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_rotor_dont_mangle_8(1.0);
		B = random_trivector_dont_mangle_7(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using rc()
		gmvC = new mv(rc(A, B));
		
		// perform GMV version 
		gmvD = new mv(rc(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("rc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_hip_dont_mangle_554(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	trivector B;
//	vector C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_rotor_dont_mangle_8(1.0);
		B = random_trivector_dont_mangle_7(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using hip()
		gmvC = new mv(hip(A, B));
		
		// perform GMV version 
		gmvD = new mv(hip(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("hip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_mhip_dont_mangle_555(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	trivector B;
//	oddVersor C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_rotor_dont_mangle_8(1.0);
		B = random_trivector_dont_mangle_7(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using mhip()
		gmvC = new mv(mhip(A, B));
		
		// perform GMV version 
		gmvD = new mv(mhip(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("mhip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_sp_dont_mangle_556(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	e1_t A;
	I3_t B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_e1_t_dont_mangle_10(1.0);
		B = random_I3_t_dont_mangle_13(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using sp()
		gmvC = new mv(sp(A, B));
		
		// perform GMV version 
		gmvD = new mv(sp(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("sp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_lc_dont_mangle_557(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	I3_t A;
	e3_t B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_I3_t_dont_mangle_13(1.0);
		B = random_e3_t_dont_mangle_75(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using lc()
		gmvC = new mv(lc(A, B));
		
		// perform GMV version 
		gmvD = new mv(lc(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("lc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_rc_dont_mangle_558(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	e1_t A;
	e1_t B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_e1_t_dont_mangle_10(1.0);
		B = random_e1_t_dont_mangle_10(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using rc()
		gmvC = new mv(rc(A, B));
		
		// perform GMV version 
		gmvD = new mv(rc(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("rc() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_hip_dont_mangle_559(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	e2_t A;
	I3_t B;
//	bivector C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_e2_t_dont_mangle_11(1.0);
		B = random_I3_t_dont_mangle_13(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using hip()
		gmvC = new mv(hip(A, B));
		
		// perform GMV version 
		gmvD = new mv(hip(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("hip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_mhip_dont_mangle_560(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	I3_t A;
	I3_t B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_I3_t_dont_mangle_13(1.0);
		B = random_I3_t_dont_mangle_13(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using mhip()
		gmvC = new mv(mhip(A, B));
		
		// perform GMV version 
		gmvD = new mv(mhip(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("mhip() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_norm_dont_mangle_561(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, reverseA, tmp;
	
	int i, g;
	int basisVectorBitmap = -1;
	double s;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		s = genrand();
		g = (int)(genrand() * 3.5);
		A = random_blade_dont_mangle_23_returns_mv(s, g, basisVectorBitmap);
		
		// compute norm
		n1 = norm(A);
		
		// compute norm manually (A . reverse(A))
		reverseA = reverse(A);
		tmp = gp(A, reverseA);
		n2 = tmp.get_scalar();
		n2 = (double)Math.Sqrt(Math.Abs(n2));
		
		// check if equal to original:
		if (Math.Abs(n1 - n2) > 1E-14) {
			Console.WriteLine("norm() test failed (difference = " + (double)Math.Abs(n1 - n2) + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_norm_dont_mangle_564(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	mv gmvA, reverseA, tmp;
	
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_vector_dont_mangle_2(genrand());
		
		gmvA = new mv(A);
		
		n1 = norm(A);
		
		// compute norm manually (A . reverse(A))
		reverseA = reverse(gmvA);
		tmp = gp(gmvA, reverseA);
		n2 = tmp.get_scalar();
		n2 = (double)Math.Sqrt(Math.Abs(n2));
		
		// check if equal to original:
		if (Math.Abs(n1 - n2) > 1E-13) {
			Console.WriteLine("norm() test failed (difference = " + (double)Math.Abs(n1 - n2) + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_norm_dont_mangle_566(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	mv gmvA, reverseA, tmp;
	
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_bivector_dont_mangle_4(genrand());
		
		gmvA = new mv(A);
		
		n1 = norm(A);
		
		// compute norm manually (A . reverse(A))
		reverseA = reverse(gmvA);
		tmp = gp(gmvA, reverseA);
		n2 = tmp.get_scalar();
		n2 = (double)Math.Sqrt(Math.Abs(n2));
		
		// check if equal to original:
		if (Math.Abs(n1 - n2) > 1E-13) {
			Console.WriteLine("norm() test failed (difference = " + (double)Math.Abs(n1 - n2) + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_norm_dont_mangle_567(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	trivector A;
	mv gmvA, reverseA, tmp;
	
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_trivector_dont_mangle_7(genrand());
		
		gmvA = new mv(A);
		
		n1 = norm(A);
		
		// compute norm manually (A . reverse(A))
		reverseA = reverse(gmvA);
		tmp = gp(gmvA, reverseA);
		n2 = tmp.get_scalar();
		n2 = (double)Math.Sqrt(Math.Abs(n2));
		
		// check if equal to original:
		if (Math.Abs(n1 - n2) > 1E-13) {
			Console.WriteLine("norm() test failed (difference = " + (double)Math.Abs(n1 - n2) + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_norm_dont_mangle_563(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	mv gmvA, reverseA, tmp;
	
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_rotor_dont_mangle_8(genrand());
		
		gmvA = new mv(A);
		
		n1 = norm(A);
		
		// compute norm manually (A . reverse(A))
		reverseA = reverse(gmvA);
		tmp = gp(gmvA, reverseA);
		n2 = tmp.get_scalar();
		n2 = (double)Math.Sqrt(Math.Abs(n2));
		
		// check if equal to original:
		if (Math.Abs(n1 - n2) > 1E-13) {
			Console.WriteLine("norm() test failed (difference = " + (double)Math.Abs(n1 - n2) + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_norm_dont_mangle_565(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	e1_t A;
	mv gmvA, reverseA, tmp;
	
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_e1_t_dont_mangle_10(genrand());
		
		gmvA = new mv(A);
		
		n1 = norm(A);
		
		// compute norm manually (A . reverse(A))
		reverseA = reverse(gmvA);
		tmp = gp(gmvA, reverseA);
		n2 = tmp.get_scalar();
		n2 = (double)Math.Sqrt(Math.Abs(n2));
		
		// check if equal to original:
		if (Math.Abs(n1 - n2) > 1E-13) {
			Console.WriteLine("norm() test failed (difference = " + (double)Math.Abs(n1 - n2) + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_norm_dont_mangle_562(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	e3_t A;
	mv gmvA, reverseA, tmp;
	
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_e3_t_dont_mangle_75(genrand());
		
		gmvA = new mv(A);
		
		n1 = norm(A);
		
		// compute norm manually (A . reverse(A))
		reverseA = reverse(gmvA);
		tmp = gp(gmvA, reverseA);
		n2 = tmp.get_scalar();
		n2 = (double)Math.Sqrt(Math.Abs(n2));
		
		// check if equal to original:
		if (Math.Abs(n1 - n2) > 1E-13) {
			Console.WriteLine("norm() test failed (difference = " + (double)Math.Abs(n1 - n2) + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_norm_dont_mangle_568(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	I3_t A;
	mv gmvA, reverseA, tmp;
	
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_I3_t_dont_mangle_13(genrand());
		
		gmvA = new mv(A);
		
		n1 = norm(A);
		
		// compute norm manually (A . reverse(A))
		reverseA = reverse(gmvA);
		tmp = gp(gmvA, reverseA);
		n2 = tmp.get_scalar();
		n2 = (double)Math.Sqrt(Math.Abs(n2));
		
		// check if equal to original:
		if (Math.Abs(n1 - n2) > 1E-13) {
			Console.WriteLine("norm() test failed (difference = " + (double)Math.Abs(n1 - n2) + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_norm2_dont_mangle_569(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, reverseA, tmp;
	
	int i, g;
	int basisVectorBitmap = -1;
	double s;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		s = genrand();
		g = (int)(genrand() * 3.5);
		A = random_blade_dont_mangle_23_returns_mv(s, g, basisVectorBitmap);
		
		// compute norm
		n1 = norm2(A);
		
		// compute norm manually (A . reverse(A))
		reverseA = reverse(A);
		tmp = gp(A, reverseA);
		n2 = tmp.get_scalar();
		
		// check if equal to original:
		if (Math.Abs(n1 - n2) > 1E-14) {
			Console.WriteLine("norm2() test failed (difference = " + (double)Math.Abs(n1 - n2) + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_norm2_dont_mangle_570(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	mv gmvA, reverseA, tmp;
	
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_vector_dont_mangle_2(genrand());
		
		gmvA = new mv(A);
		
		n1 = norm2(A);
		
		// compute norm manually (A . reverse(A))
		reverseA = reverse(gmvA);
		tmp = gp(gmvA, reverseA);
		n2 = tmp.get_scalar();
		
		// check if equal to original:
		if (Math.Abs(n1 - n2) > 1E-13) {
			Console.WriteLine("norm2() test failed (difference = " + (double)Math.Abs(n1 - n2) + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_norm2_dont_mangle_571(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	mv gmvA, reverseA, tmp;
	
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_bivector_dont_mangle_4(genrand());
		
		gmvA = new mv(A);
		
		n1 = norm2(A);
		
		// compute norm manually (A . reverse(A))
		reverseA = reverse(gmvA);
		tmp = gp(gmvA, reverseA);
		n2 = tmp.get_scalar();
		
		// check if equal to original:
		if (Math.Abs(n1 - n2) > 1E-13) {
			Console.WriteLine("norm2() test failed (difference = " + (double)Math.Abs(n1 - n2) + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_norm2_dont_mangle_572(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	trivector A;
	mv gmvA, reverseA, tmp;
	
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_trivector_dont_mangle_7(genrand());
		
		gmvA = new mv(A);
		
		n1 = norm2(A);
		
		// compute norm manually (A . reverse(A))
		reverseA = reverse(gmvA);
		tmp = gp(gmvA, reverseA);
		n2 = tmp.get_scalar();
		
		// check if equal to original:
		if (Math.Abs(n1 - n2) > 1E-13) {
			Console.WriteLine("norm2() test failed (difference = " + (double)Math.Abs(n1 - n2) + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_norm2_dont_mangle_573(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	mv gmvA, reverseA, tmp;
	
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_rotor_dont_mangle_8(genrand());
		
		gmvA = new mv(A);
		
		n1 = norm2(A);
		
		// compute norm manually (A . reverse(A))
		reverseA = reverse(gmvA);
		tmp = gp(gmvA, reverseA);
		n2 = tmp.get_scalar();
		
		// check if equal to original:
		if (Math.Abs(n1 - n2) > 1E-13) {
			Console.WriteLine("norm2() test failed (difference = " + (double)Math.Abs(n1 - n2) + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_norm2_dont_mangle_574(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	e1_t A;
	mv gmvA, reverseA, tmp;
	
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_e1_t_dont_mangle_10(genrand());
		
		gmvA = new mv(A);
		
		n1 = norm2(A);
		
		// compute norm manually (A . reverse(A))
		reverseA = reverse(gmvA);
		tmp = gp(gmvA, reverseA);
		n2 = tmp.get_scalar();
		
		// check if equal to original:
		if (Math.Abs(n1 - n2) > 1E-13) {
			Console.WriteLine("norm2() test failed (difference = " + (double)Math.Abs(n1 - n2) + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_norm2_dont_mangle_575(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	e3_t A;
	mv gmvA, reverseA, tmp;
	
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_e3_t_dont_mangle_75(genrand());
		
		gmvA = new mv(A);
		
		n1 = norm2(A);
		
		// compute norm manually (A . reverse(A))
		reverseA = reverse(gmvA);
		tmp = gp(gmvA, reverseA);
		n2 = tmp.get_scalar();
		
		// check if equal to original:
		if (Math.Abs(n1 - n2) > 1E-13) {
			Console.WriteLine("norm2() test failed (difference = " + (double)Math.Abs(n1 - n2) + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_norm2_dont_mangle_576(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	I3_t A;
	mv gmvA, reverseA, tmp;
	
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_I3_t_dont_mangle_13(genrand());
		
		gmvA = new mv(A);
		
		n1 = norm2(A);
		
		// compute norm manually (A . reverse(A))
		reverseA = reverse(gmvA);
		tmp = gp(gmvA, reverseA);
		n2 = tmp.get_scalar();
		
		// check if equal to original:
		if (Math.Abs(n1 - n2) > 1E-13) {
			Console.WriteLine("norm2() test failed (difference = " + (double)Math.Abs(n1 - n2) + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_op_dont_mangle_577_1(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, C, D, E, dif;
	int i, ga, gb, gd;
	double s;
	int basisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		ga = (int)(genrand() * 3.5);
		A = random_blade_dont_mangle_23_returns_mv(s, ga, basisVectorBitmap);
		
		s = genrand();
		gb = (int)(genrand() * 3.5);
		B = random_blade_dont_mangle_23_returns_mv(s, gb, basisVectorBitmap);
		
		// compute product using op()
		C = new mv(op(A, B));
		
		// simulate product using geometric product & grade part selection
		D = gp(A, B);
		gd = (ga > gb) ? ga - gb : gb - ga;
		E = extractGrade(D,  Grades[ga + gb]);

		// check if equal:
		dif = subtract(C, E);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("op() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_op_dont_mangle_577_2(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B;
	int i, ga;
	double s;
	int basisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		ga = (int)(genrand() * 3.5);
		if (ga == 0) continue; // do not perform this test for grade 0
		A = random_blade_dont_mangle_23_returns_mv(s, ga, basisVectorBitmap);
		
		// compute A ^ A (should be zero)
		B = op(A, A);
		
		// check if B is zero:
		if (B.LargestCoordinate() > 1E-13) {
			Console.WriteLine("op() test failed (largestCoordinate = " + (double)B.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_op_dont_mangle_578(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	vector B;
//	bivector C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_vector_dont_mangle_2(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using op()
		gmvC = new mv(op(A, B));
		
		// perform GMV version 
		gmvD = new mv(op(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("op() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_op_dont_mangle_579(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	bivector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_bivector_dont_mangle_4(1.0);
		B = random_bivector_dont_mangle_4(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using op()
		gmvC = new mv(op(A, B));
		
		// perform GMV version 
		gmvD = new mv(op(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("op() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_op_dont_mangle_580(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	rotor B;
//	oddVersor C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_rotor_dont_mangle_8(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using op()
		gmvC = new mv(op(A, B));
		
		// perform GMV version 
		gmvD = new mv(op(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("op() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_op_dont_mangle_581(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	trivector B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(1.0);
		B = random_trivector_dont_mangle_7(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using op()
		gmvC = new mv(op(A, B));
		
		// perform GMV version 
		gmvD = new mv(op(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("op() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_op_dont_mangle_582(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	e1_t A;
	e1_t B;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_e1_t_dont_mangle_10(1.0);
		B = random_e1_t_dont_mangle_10(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using op()
		gmvC = new mv(op(A, B));
		
		// perform GMV version 
		gmvD = new mv(op(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("op() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_op_dont_mangle_583(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	e1_t A;
	e2_t B;
//	bivector C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_e1_t_dont_mangle_10(1.0);
		B = random_e2_t_dont_mangle_11(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using op()
		gmvC = new mv(op(A, B));
		
		// perform GMV version 
		gmvD = new mv(op(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("op() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_op_dont_mangle_584(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	e2_t A;
	e3_t B;
//	bivector C;
	

	mv gmvA, gmvB, gmvC, gmvD, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_e2_t_dont_mangle_11(1.0);
		B = random_e3_t_dont_mangle_75(1.0);
		
		// convert smv to GMV
		gmvA = new mv(A);
		gmvB = new mv(B);
				
		// compute product using op()
		gmvC = new mv(op(A, B));
		
		// perform GMV version 
		gmvD = new mv(op(gmvA, gmvB));

		// check if equal:
		dif = subtract(gmvC, gmvD);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("op() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_exp_dont_mangle_586(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, expA, sinhA, coshA, S, dif, tmp1; //, tmp2;
	int i, g;
	int basisVectorBitmap = -1;
	double s;
	int order = 12;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade of grade 2
		s = 2.0 * genrand();
		g = 2;
		A = random_blade_dont_mangle_23_returns_mv(s, g, basisVectorBitmap);
		
		if (genrand() > 0.5) { // make sure that 'A' is not always a blade
			s = genrand();
			tmp1 = random_blade_dont_mangle_23_returns_mv(s, g, basisVectorBitmap);	
			A = add(A, tmp1);
			//A = tmp2;
		}

		// apply sinh, cosh, exp functions
		expA = exp(A, order);
		sinhA = sinh(A, order);
		coshA = cosh(A, order);
		
		// sum sinh and cosh
		S = add(coshA, sinhA);
		
		// test that sinh+cosh == exp
		dif = subtract(expA, S);
		if (dif.LargestCoordinate() > 0.00031622776601683783) { // sinh, cosh precision is very low
			Console.WriteLine("exp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_exp_dont_mangle_587(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	rotor B;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	int order = 12;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_bivector_dont_mangle_4(genrand());
		
		// compute sin, cos or exp
		B = exp(A);
		gmvB = new mv(B);
		
		// compute sin, cos or exp using GMV
		gmvA = new mv(A); 
		gmvC = exp(gmvA, order);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("exp() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_cosh_dont_mangle_589(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, expA, sinhA, coshA, S, dif, tmp1; //, tmp2;
	int i, g;
	int basisVectorBitmap = -1;
	double s;
	int order = 12;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade of grade 2
		s = 2.0 * genrand();
		g = 2;
		A = random_blade_dont_mangle_23_returns_mv(s, g, basisVectorBitmap);
		
		if (genrand() > 0.5) { // make sure that 'A' is not always a blade
			s = genrand();
			tmp1 = random_blade_dont_mangle_23_returns_mv(s, g, basisVectorBitmap);	
			A = add(A, tmp1);
			//A = tmp2;
		}

		// apply sinh, cosh, exp functions
		expA = exp(A, order);
		sinhA = sinh(A, order);
		coshA = cosh(A, order);
		
		// sum sinh and cosh
		S = add(coshA, sinhA);
		
		// test that sinh+cosh == exp
		dif = subtract(expA, S);
		if (dif.LargestCoordinate() > 0.00031622776601683783) { // sinh, cosh precision is very low
			Console.WriteLine("cosh() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_cosh_dont_mangle_588(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	double B;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	int order = 12;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_bivector_dont_mangle_4(genrand());
		
		// compute sin, cos or exp
		B = cosh(A);
		gmvB = new mv(B);
		
		// compute sin, cos or exp using GMV
		gmvA = new mv(A); 
		gmvC = cosh(gmvA, order);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("cosh() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_sinh_dont_mangle_590(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, expA, sinhA, coshA, S, dif, tmp1; //, tmp2;
	int i, g;
	int basisVectorBitmap = -1;
	double s;
	int order = 12;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade of grade 2
		s = 2.0 * genrand();
		g = 2;
		A = random_blade_dont_mangle_23_returns_mv(s, g, basisVectorBitmap);
		
		if (genrand() > 0.5) { // make sure that 'A' is not always a blade
			s = genrand();
			tmp1 = random_blade_dont_mangle_23_returns_mv(s, g, basisVectorBitmap);	
			A = add(A, tmp1);
			//A = tmp2;
		}

		// apply sinh, cosh, exp functions
		expA = exp(A, order);
		sinhA = sinh(A, order);
		coshA = cosh(A, order);
		
		// sum sinh and cosh
		S = add(coshA, sinhA);
		
		// test that sinh+cosh == exp
		dif = subtract(expA, S);
		if (dif.LargestCoordinate() > 0.00031622776601683783) { // sinh, cosh precision is very low
			Console.WriteLine("sinh() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_sinh_dont_mangle_591(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	bivector B;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	int order = 12;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_bivector_dont_mangle_4(genrand());
		
		// compute sin, cos or exp
		B = sinh(A);
		gmvB = new mv(B);
		
		// compute sin, cos or exp using GMV
		gmvA = new mv(A); 
		gmvC = sinh(gmvA, order);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("sinh() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_cos_dont_mangle_593(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	double B;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	int order = 12;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_bivector_dont_mangle_4(genrand());
		
		// compute sin, cos or exp
		B = cos(A);
		gmvB = new mv(B);
		
		// compute sin, cos or exp using GMV
		gmvA = new mv(A); 
		gmvC = cos(gmvA, order);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("cos() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_sin_dont_mangle_595(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	bivector B;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	int order = 12;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_bivector_dont_mangle_4(genrand());
		
		// compute sin, cos or exp
		B = sin(A);
		gmvB = new mv(B);
		
		// compute sin, cos or exp using GMV
		gmvA = new mv(A); 
		gmvC = sin(gmvA, order);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-12) {
			Console.WriteLine("sin() test failed (largestCoordinate = " + (double)dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_negate_dont_mangle_596(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, GA, GB;
	int i, c, n, g;
	int basisVectorBitmap = -1;
	double s, dif;
	int[] signTable = new int[]{-1, -1, -1, -1};
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		s = genrand();
		g = (int)(genrand() * 3.5);
		A = random_versor_dont_mangle_1_returns_mv(s, g, basisVectorBitmap);

		// apply function
		B = negate(A);
		
		// check grade
		for (n = 0; n <= 3; n++) {
			GA = extractGrade(A, Grades[n]);
			GB = extractGrade(B, Grades[n]);
			
			// check if grade usage matches
			if (GA.gu() != GB.gu()) {
				Console.WriteLine("negate() test failed (grade usage does not match)");
				return 0; // failure
			}
			
			// check each coordinate of each groups which is still present after the grade selection
			for (int m = 0; m < 4; m++) {
				if (GA.m_c[m] != null) {
					for (c = 0; c < GA.m_c[m].Length; c++) {
						dif = (double)Math.Abs(GA.m_c[m][c] * (double)signTable[n] - GB.m_c[m][c]);
						if (dif > 1E-14) {
							Console.WriteLine("negate() test failed (dif = " + dif + ")");
							return 0; // failure
						}
					}
				}
			}
		}
		
	}
	return 1; // success
}

static int test_cliffordConjugate_dont_mangle_597(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, GA, GB;
	int i, c, n, g;
	int basisVectorBitmap = -1;
	double s, dif;
	int[] signTable = new int[]{1, -1, -1, 1};
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		s = genrand();
		g = (int)(genrand() * 3.5);
		A = random_versor_dont_mangle_1_returns_mv(s, g, basisVectorBitmap);

		// apply function
		B = cliffordConjugate(A);
		
		// check grade
		for (n = 0; n <= 3; n++) {
			GA = extractGrade(A, Grades[n]);
			GB = extractGrade(B, Grades[n]);
			
			// check if grade usage matches
			if (GA.gu() != GB.gu()) {
				Console.WriteLine("cliffordConjugate() test failed (grade usage does not match)");
				return 0; // failure
			}
			
			// check each coordinate of each groups which is still present after the grade selection
			for (int m = 0; m < 4; m++) {
				if (GA.m_c[m] != null) {
					for (c = 0; c < GA.m_c[m].Length; c++) {
						dif = (double)Math.Abs(GA.m_c[m][c] * (double)signTable[n] - GB.m_c[m][c]);
						if (dif > 1E-14) {
							Console.WriteLine("cliffordConjugate() test failed (dif = " + dif + ")");
							return 0; // failure
						}
					}
				}
			}
		}
		
	}
	return 1; // success
}

static int test_gradeInvolution_dont_mangle_598(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, GA, GB;
	int i, c, n, g;
	int basisVectorBitmap = -1;
	double s, dif;
	int[] signTable = new int[]{1, -1, 1, -1};
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		s = genrand();
		g = (int)(genrand() * 3.5);
		A = random_versor_dont_mangle_1_returns_mv(s, g, basisVectorBitmap);

		// apply function
		B = gradeInvolution(A);
		
		// check grade
		for (n = 0; n <= 3; n++) {
			GA = extractGrade(A, Grades[n]);
			GB = extractGrade(B, Grades[n]);
			
			// check if grade usage matches
			if (GA.gu() != GB.gu()) {
				Console.WriteLine("gradeInvolution() test failed (grade usage does not match)");
				return 0; // failure
			}
			
			// check each coordinate of each groups which is still present after the grade selection
			for (int m = 0; m < 4; m++) {
				if (GA.m_c[m] != null) {
					for (c = 0; c < GA.m_c[m].Length; c++) {
						dif = (double)Math.Abs(GA.m_c[m][c] * (double)signTable[n] - GB.m_c[m][c]);
						if (dif > 1E-14) {
							Console.WriteLine("gradeInvolution() test failed (dif = " + dif + ")");
							return 0; // failure
						}
					}
				}
			}
		}
		
	}
	return 1; // success
}

static int test_reverse_dont_mangle_599(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, B, GA, GB;
	int i, c, n, g;
	int basisVectorBitmap = -1;
	double s, dif;
	int[] signTable = new int[]{1, 1, -1, -1};
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		s = genrand();
		g = (int)(genrand() * 3.5);
		A = random_versor_dont_mangle_1_returns_mv(s, g, basisVectorBitmap);

		// apply function
		B = reverse(A);
		
		// check grade
		for (n = 0; n <= 3; n++) {
			GA = extractGrade(A, Grades[n]);
			GB = extractGrade(B, Grades[n]);
			
			// check if grade usage matches
			if (GA.gu() != GB.gu()) {
				Console.WriteLine("reverse() test failed (grade usage does not match)");
				return 0; // failure
			}
			
			// check each coordinate of each groups which is still present after the grade selection
			for (int m = 0; m < 4; m++) {
				if (GA.m_c[m] != null) {
					for (c = 0; c < GA.m_c[m].Length; c++) {
						dif = (double)Math.Abs(GA.m_c[m][c] * (double)signTable[n] - GB.m_c[m][c]);
						if (dif > 1E-14) {
							Console.WriteLine("reverse() test failed (dif = " + dif + ")");
							return 0; // failure
						}
					}
				}
			}
		}
		
	}
	return 1; // success
}

static int test_negate_dont_mangle_600(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_vector_dont_mangle_2(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(negate(A));
		
		// compute via GMV
		gmvC = negate(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("negate() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_cliffordConjugate_dont_mangle_601(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_vector_dont_mangle_2(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(cliffordConjugate(A));
		
		// compute via GMV
		gmvC = cliffordConjugate(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("cliffordConjugate() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_gradeInvolution_dont_mangle_602(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_vector_dont_mangle_2(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(gradeInvolution(A));
		
		// compute via GMV
		gmvC = gradeInvolution(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("gradeInvolution() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_reverse_dont_mangle_603(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_vector_dont_mangle_2(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(reverse(A));
		
		// compute via GMV
		gmvC = reverse(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("reverse() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_negate_dont_mangle_604(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_bivector_dont_mangle_4(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(negate(A));
		
		// compute via GMV
		gmvC = negate(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("negate() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_cliffordConjugate_dont_mangle_605(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_bivector_dont_mangle_4(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(cliffordConjugate(A));
		
		// compute via GMV
		gmvC = cliffordConjugate(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("cliffordConjugate() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_gradeInvolution_dont_mangle_606(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_bivector_dont_mangle_4(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(gradeInvolution(A));
		
		// compute via GMV
		gmvC = gradeInvolution(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("gradeInvolution() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_reverse_dont_mangle_607(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_bivector_dont_mangle_4(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(reverse(A));
		
		// compute via GMV
		gmvC = reverse(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("reverse() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_negate_dont_mangle_608(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	trivector A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_trivector_dont_mangle_7(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(negate(A));
		
		// compute via GMV
		gmvC = negate(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("negate() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_cliffordConjugate_dont_mangle_609(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	trivector A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_trivector_dont_mangle_7(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(cliffordConjugate(A));
		
		// compute via GMV
		gmvC = cliffordConjugate(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("cliffordConjugate() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_gradeInvolution_dont_mangle_610(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	trivector A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_trivector_dont_mangle_7(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(gradeInvolution(A));
		
		// compute via GMV
		gmvC = gradeInvolution(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("gradeInvolution() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_reverse_dont_mangle_611(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	trivector A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_trivector_dont_mangle_7(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(reverse(A));
		
		// compute via GMV
		gmvC = reverse(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("reverse() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_negate_dont_mangle_612(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_rotor_dont_mangle_8(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(negate(A));
		
		// compute via GMV
		gmvC = negate(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("negate() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_cliffordConjugate_dont_mangle_613(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_rotor_dont_mangle_8(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(cliffordConjugate(A));
		
		// compute via GMV
		gmvC = cliffordConjugate(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("cliffordConjugate() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_gradeInvolution_dont_mangle_614(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_rotor_dont_mangle_8(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(gradeInvolution(A));
		
		// compute via GMV
		gmvC = gradeInvolution(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("gradeInvolution() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_reverse_dont_mangle_615(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_rotor_dont_mangle_8(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(reverse(A));
		
		// compute via GMV
		gmvC = reverse(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("reverse() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_negate_dont_mangle_616(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	e1_t A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_e1_t_dont_mangle_10(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(negate(A));
		
		// compute via GMV
		gmvC = negate(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("negate() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_cliffordConjugate_dont_mangle_617(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	e2_t A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_e2_t_dont_mangle_11(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(cliffordConjugate(A));
		
		// compute via GMV
		gmvC = cliffordConjugate(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("cliffordConjugate() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_gradeInvolution_dont_mangle_618(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	e3_t A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_e3_t_dont_mangle_75(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(gradeInvolution(A));
		
		// compute via GMV
		gmvC = gradeInvolution(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("gradeInvolution() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_reverse_dont_mangle_619(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	I3_t A;
	mv gmvA, gmvB, gmvC, dif;
	int i;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random Smv
		A = random_I3_t_dont_mangle_13(genrand());
		gmvA = new mv(A);
		
		gmvB = new mv(reverse(A));
		
		// compute via GMV
		gmvC = reverse(gmvA);
		
		// check if equal:
		dif = subtract(gmvC, gmvB);
		if (dif.LargestCoordinate() > 1E-14) {
			Console.WriteLine("reverse() test failed (largestCoordinate = " + dif.LargestCoordinate() + ")");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_unit_dont_mangle_624(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A, UA, RUA;
	int i;
	int basisVectorBitmap = -1;
	double n;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		A = random_blade_dont_mangle_23_returns_mv(genrand(), (int)(genrand() * 3.5), basisVectorBitmap);
		
		UA = unit(A);
		RUA = reverse(UA);
		n = sp(RUA, UA);
		
		if ((double)(Math.Abs(n) - 1.0) > 1E-12) {
			Console.WriteLine("unit() test failed (|norm|-1 = " + (double)(Math.Abs((double)n) - 1.0) + ")");
			return 0; // failure
		}

	}
	return 1; // success
}

static int test_unit_dont_mangle_626(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	vector A;
	vector UA;
	mv gmvUA, RUA;
	int i;
	double n;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_vector_dont_mangle_2(genrand());
		
		UA = unit(A);
		
		gmvUA = new mv(UA);
		
		RUA = reverse(gmvUA);
		n = sp(RUA, gmvUA);
		
		if ((double)(Math.Abs(n) - 1.0) > 1E-12) {
			Console.WriteLine("unit() test failed (|norm|-1 = " + (double)(Math.Abs((double)n) - 1.0) + ")");
			return 0; // failure
		}

	}
	return 1; // success
}

static int test_unit_dont_mangle_628(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	bivector A;
	bivector UA;
	mv gmvUA, RUA;
	int i;
	double n;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_bivector_dont_mangle_4(genrand());
		
		UA = unit(A);
		
		gmvUA = new mv(UA);
		
		RUA = reverse(gmvUA);
		n = sp(RUA, gmvUA);
		
		if ((double)(Math.Abs(n) - 1.0) > 1E-12) {
			Console.WriteLine("unit() test failed (|norm|-1 = " + (double)(Math.Abs((double)n) - 1.0) + ")");
			return 0; // failure
		}

	}
	return 1; // success
}

static int test_unit_dont_mangle_627(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	trivector A;
	trivector UA;
	mv gmvUA, RUA;
	int i;
	double n;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_trivector_dont_mangle_7(genrand());
		
		UA = unit(A);
		
		gmvUA = new mv(UA);
		
		RUA = reverse(gmvUA);
		n = sp(RUA, gmvUA);
		
		if ((double)(Math.Abs(n) - 1.0) > 1E-12) {
			Console.WriteLine("unit() test failed (|norm|-1 = " + (double)(Math.Abs((double)n) - 1.0) + ")");
			return 0; // failure
		}

	}
	return 1; // success
}

static int test_unit_dont_mangle_625(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	rotor A;
	rotor UA;
	mv gmvUA, RUA;
	int i;
	double n;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_rotor_dont_mangle_8(genrand());
		
		UA = unit(A);
		
		gmvUA = new mv(UA);
		
		RUA = reverse(gmvUA);
		n = sp(RUA, gmvUA);
		
		if ((double)(Math.Abs(n) - 1.0) > 1E-12) {
			Console.WriteLine("unit() test failed (|norm|-1 = " + (double)(Math.Abs((double)n) - 1.0) + ")");
			return 0; // failure
		}

	}
	return 1; // success
}

static int test_unit_dont_mangle_629(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	oddVersor A;
	oddVersor UA;
	mv gmvUA, RUA;
	int i;
	double n;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_oddVersor_dont_mangle_93(genrand());
		
		UA = unit(A);
		
		gmvUA = new mv(UA);
		
		RUA = reverse(gmvUA);
		n = sp(RUA, gmvUA);
		
		if ((double)(Math.Abs(n) - 1.0) > 1E-12) {
			Console.WriteLine("unit() test failed (|norm|-1 = " + (double)(Math.Abs((double)n) - 1.0) + ")");
			return 0; // failure
		}

	}
	return 1; // success
}

static int test_unit_dont_mangle_630(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	e1_t A;
	e1_t UA;
	mv gmvUA, RUA;
	int i;
	double n;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_e1_t_dont_mangle_10(genrand());
		
		UA = unit(A);
		
		gmvUA = new mv(UA);
		
		RUA = reverse(gmvUA);
		n = sp(RUA, gmvUA);
		
		if ((double)(Math.Abs(n) - 1.0) > 1E-12) {
			Console.WriteLine("unit() test failed (|norm|-1 = " + (double)(Math.Abs((double)n) - 1.0) + ")");
			return 0; // failure
		}

	}
	return 1; // success
}

static int test_unit_dont_mangle_631(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	e2_t A;
	e2_t UA;
	mv gmvUA, RUA;
	int i;
	double n;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_e2_t_dont_mangle_11(genrand());
		
		UA = unit(A);
		
		gmvUA = new mv(UA);
		
		RUA = reverse(gmvUA);
		n = sp(RUA, gmvUA);
		
		if ((double)(Math.Abs(n) - 1.0) > 1E-12) {
			Console.WriteLine("unit() test failed (|norm|-1 = " + (double)(Math.Abs((double)n) - 1.0) + ")");
			return 0; // failure
		}

	}
	return 1; // success
}

static int test_unit_dont_mangle_632(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	I3_t A;
	I3_t UA;
	mv gmvUA, RUA;
	int i;
	double n;
	
	for (i = 0; i < NB_LOOPS; i++) {
		A = random_I3_t_dont_mangle_13(genrand());
		
		UA = unit(A);
		
		gmvUA = new mv(UA);
		
		RUA = reverse(gmvUA);
		n = sp(RUA, gmvUA);
		
		if ((double)(Math.Abs(n) - 1.0) > 1E-12) {
			Console.WriteLine("unit() test failed (|norm|-1 = " + (double)(Math.Abs((double)n) - 1.0) + ")");
			return 0; // failure
		}

	}
	return 1; // success
}

static int test_versorInverse_dont_mangle_633(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv V, VI, VVI, VIV, X, Y;
	int i;
	int basisVectorBitmap = -1;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random versor
		V = random_versor_dont_mangle_1_returns_mv(genrand() + 0.5, (int)(genrand() * 3.5), basisVectorBitmap);
		
		// avoid 'near'-singular versors
		if (V.LargestCoordinate() > 2.0)
			continue;
		
		// compute inverse
		VI = versorInverse(V);
		
		// compute (inverse(V) V) and (V inverse(V))
		VIV = gp(VI, V);
		VVI = gp(V, VI);
		
		// check that scalar parts are 1
		n1 = VIV.get_scalar();
		n2 = VVI.get_scalar();
		
		if (Math.Abs(n1 - 1.0) > 1E-11) {
			Console.WriteLine("versorInverse() test failed |VI . V - 1|= " + (double)Math.Abs(n1 - 1.0) + ")");
			return 0; // failure
		}
		
		if (Math.Abs(n2 - 1.0) > 1E-11) {
			Console.WriteLine("versorInverse() test failed ( |V . VI - 1| = " + (double)Math.Abs(n2 - 1.0) + ")");
			return 0; // failure
		}
		
		// check that other grade parts are zero:
		X = extractGrade(VIV, GroupBitmap.ALL_GROUPS ^ GroupBitmap.GROUP_0);
		Y = extractGrade(VVI, GroupBitmap.ALL_GROUPS ^ GroupBitmap.GROUP_0);
		
		if (X.LargestCoordinate() > 1E-11) {
			Console.WriteLine("versorInverse() test failed (largestCoordinate of VIV = " + (double)X.LargestCoordinate() + ")");
			return 0; // failure
		}
		
		if (Y.LargestCoordinate() > 1E-11) {
			Console.WriteLine("versorInverse() test failed (largestCoordinate of VVI = " + (double)Y.LargestCoordinate() + ")");
			return 0; // failure
		}
		
	}
	return 1; // success
}

static int test_versorInverse_dont_mangle_634(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	vector V;
	vector VI;
	mv gmvV, gmvVI, VVI, VIV, X, Y;
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		V = random_vector_dont_mangle_2(genrand() + 0.5);
		
		// avoid 'near'-singular versors
		if (V.LargestCoordinate() > 2.0)
			continue;
		
		// compute inverse
		VI = versorInverse(V);
		
		// convert to GMV
		gmvV = new mv(V);
		gmvVI = new mv(VI);
		
		// compute (inverse(V) V) and (V inverse(V))
		VIV = gp(gmvVI, gmvV);
		VVI = gp(gmvV, gmvVI);
		
		// check that scalar parts are 1
		n1 = VIV.get_scalar();
		n2 = VVI.get_scalar();
		
		if (Math.Abs(n1 - 1.0) > 1E-11) {
			Console.WriteLine("versorInverse() test failed |VI . V - 1|= " + (double)Math.Abs(n1 - 1.0) + ")");
			return 0; // failure
		}
		
		if (Math.Abs(n2 - 1.0) > 1E-11) {
			Console.WriteLine("versorInverse() test failed ( |V . VI - 1| = " + (double)Math.Abs(n2 - 1.0) + ")");
			return 0; // failure
		}
		
		// check that other grade parts are zero:
		X = extractGrade(VIV, GroupBitmap.ALL_GROUPS ^ GroupBitmap.GROUP_0);
		Y = extractGrade(VVI, GroupBitmap.ALL_GROUPS ^ GroupBitmap.GROUP_0);
		
		if (X.LargestCoordinate() > 1E-11) {
			Console.WriteLine("versorInverse() test failed (largestCoordinate of VIV = " + (double)X.LargestCoordinate() + ")");
			return 0; // failure
		}
		
		if (Y.LargestCoordinate() > 1E-11) {
			Console.WriteLine("versorInverse() test failed (largestCoordinate of VVI = " + (double)Y.LargestCoordinate() + ")");
			return 0; // failure
		}
		
	}
	return 1; // success
}

static int test_versorInverse_dont_mangle_635(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 3;
	bivector V;
	bivector VI;
	mv gmvV, gmvVI, VVI, VIV, X, Y;
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		V = random_bivector_dont_mangle_4(genrand() + 0.5);
		
		// avoid 'near'-singular versors
		if (V.LargestCoordinate() > 2.0)
			continue;
		
		// compute inverse
		VI = versorInverse(V);
		
		// convert to GMV
		gmvV = new mv(V);
		gmvVI = new mv(VI);
		
		// compute (inverse(V) V) and (V inverse(V))
		VIV = gp(gmvVI, gmvV);
		VVI = gp(gmvV, gmvVI);
		
		// check that scalar parts are 1
		n1 = VIV.get_scalar();
		n2 = VVI.get_scalar();
		
		if (Math.Abs(n1 - 1.0) > 1E-11) {
			Console.WriteLine("versorInverse() test failed |VI . V - 1|= " + (double)Math.Abs(n1 - 1.0) + ")");
			return 0; // failure
		}
		
		if (Math.Abs(n2 - 1.0) > 1E-11) {
			Console.WriteLine("versorInverse() test failed ( |V . VI - 1| = " + (double)Math.Abs(n2 - 1.0) + ")");
			return 0; // failure
		}
		
		// check that other grade parts are zero:
		X = extractGrade(VIV, GroupBitmap.ALL_GROUPS ^ GroupBitmap.GROUP_0);
		Y = extractGrade(VVI, GroupBitmap.ALL_GROUPS ^ GroupBitmap.GROUP_0);
		
		if (X.LargestCoordinate() > 1E-11) {
			Console.WriteLine("versorInverse() test failed (largestCoordinate of VIV = " + (double)X.LargestCoordinate() + ")");
			return 0; // failure
		}
		
		if (Y.LargestCoordinate() > 1E-11) {
			Console.WriteLine("versorInverse() test failed (largestCoordinate of VVI = " + (double)Y.LargestCoordinate() + ")");
			return 0; // failure
		}
		
	}
	return 1; // success
}

static int test_versorInverse_dont_mangle_636(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	trivector V;
	trivector VI;
	mv gmvV, gmvVI, VVI, VIV, X, Y;
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		V = random_trivector_dont_mangle_7(genrand() + 0.5);
		
		// avoid 'near'-singular versors
		if (V.LargestCoordinate() > 2.0)
			continue;
		
		// compute inverse
		VI = versorInverse(V);
		
		// convert to GMV
		gmvV = new mv(V);
		gmvVI = new mv(VI);
		
		// compute (inverse(V) V) and (V inverse(V))
		VIV = gp(gmvVI, gmvV);
		VVI = gp(gmvV, gmvVI);
		
		// check that scalar parts are 1
		n1 = VIV.get_scalar();
		n2 = VVI.get_scalar();
		
		if (Math.Abs(n1 - 1.0) > 1E-11) {
			Console.WriteLine("versorInverse() test failed |VI . V - 1|= " + (double)Math.Abs(n1 - 1.0) + ")");
			return 0; // failure
		}
		
		if (Math.Abs(n2 - 1.0) > 1E-11) {
			Console.WriteLine("versorInverse() test failed ( |V . VI - 1| = " + (double)Math.Abs(n2 - 1.0) + ")");
			return 0; // failure
		}
		
		// check that other grade parts are zero:
		X = extractGrade(VIV, GroupBitmap.ALL_GROUPS ^ GroupBitmap.GROUP_0);
		Y = extractGrade(VVI, GroupBitmap.ALL_GROUPS ^ GroupBitmap.GROUP_0);
		
		if (X.LargestCoordinate() > 1E-11) {
			Console.WriteLine("versorInverse() test failed (largestCoordinate of VIV = " + (double)X.LargestCoordinate() + ")");
			return 0; // failure
		}
		
		if (Y.LargestCoordinate() > 1E-11) {
			Console.WriteLine("versorInverse() test failed (largestCoordinate of VVI = " + (double)Y.LargestCoordinate() + ")");
			return 0; // failure
		}
		
	}
	return 1; // success
}

static int test_versorInverse_dont_mangle_637(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	rotor V;
	rotor VI;
	mv gmvV, gmvVI, VVI, VIV, X, Y;
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		V = random_rotor_dont_mangle_8(genrand() + 0.5);
		
		// avoid 'near'-singular versors
		if (V.LargestCoordinate() > 2.0)
			continue;
		
		// compute inverse
		VI = versorInverse(V);
		
		// convert to GMV
		gmvV = new mv(V);
		gmvVI = new mv(VI);
		
		// compute (inverse(V) V) and (V inverse(V))
		VIV = gp(gmvVI, gmvV);
		VVI = gp(gmvV, gmvVI);
		
		// check that scalar parts are 1
		n1 = VIV.get_scalar();
		n2 = VVI.get_scalar();
		
		if (Math.Abs(n1 - 1.0) > 1E-11) {
			Console.WriteLine("versorInverse() test failed |VI . V - 1|= " + (double)Math.Abs(n1 - 1.0) + ")");
			return 0; // failure
		}
		
		if (Math.Abs(n2 - 1.0) > 1E-11) {
			Console.WriteLine("versorInverse() test failed ( |V . VI - 1| = " + (double)Math.Abs(n2 - 1.0) + ")");
			return 0; // failure
		}
		
		// check that other grade parts are zero:
		X = extractGrade(VIV, GroupBitmap.ALL_GROUPS ^ GroupBitmap.GROUP_0);
		Y = extractGrade(VVI, GroupBitmap.ALL_GROUPS ^ GroupBitmap.GROUP_0);
		
		if (X.LargestCoordinate() > 1E-11) {
			Console.WriteLine("versorInverse() test failed (largestCoordinate of VIV = " + (double)X.LargestCoordinate() + ")");
			return 0; // failure
		}
		
		if (Y.LargestCoordinate() > 1E-11) {
			Console.WriteLine("versorInverse() test failed (largestCoordinate of VVI = " + (double)Y.LargestCoordinate() + ")");
			return 0; // failure
		}
		
	}
	return 1; // success
}

static int test_versorInverse_dont_mangle_638(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	e1_t V;
	e1_t VI;
	mv gmvV, gmvVI, VVI, VIV, X, Y;
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		V = random_e1_t_dont_mangle_10(genrand() + 0.5);
		
		// avoid 'near'-singular versors
		if (V.LargestCoordinate() > 2.0)
			continue;
		
		// compute inverse
		VI = versorInverse(V);
		
		// convert to GMV
		gmvV = new mv(V);
		gmvVI = new mv(VI);
		
		// compute (inverse(V) V) and (V inverse(V))
		VIV = gp(gmvVI, gmvV);
		VVI = gp(gmvV, gmvVI);
		
		// check that scalar parts are 1
		n1 = VIV.get_scalar();
		n2 = VVI.get_scalar();
		
		if (Math.Abs(n1 - 1.0) > 1E-11) {
			Console.WriteLine("versorInverse() test failed |VI . V - 1|= " + (double)Math.Abs(n1 - 1.0) + ")");
			return 0; // failure
		}
		
		if (Math.Abs(n2 - 1.0) > 1E-11) {
			Console.WriteLine("versorInverse() test failed ( |V . VI - 1| = " + (double)Math.Abs(n2 - 1.0) + ")");
			return 0; // failure
		}
		
		// check that other grade parts are zero:
		X = extractGrade(VIV, GroupBitmap.ALL_GROUPS ^ GroupBitmap.GROUP_0);
		Y = extractGrade(VVI, GroupBitmap.ALL_GROUPS ^ GroupBitmap.GROUP_0);
		
		if (X.LargestCoordinate() > 1E-11) {
			Console.WriteLine("versorInverse() test failed (largestCoordinate of VIV = " + (double)X.LargestCoordinate() + ")");
			return 0; // failure
		}
		
		if (Y.LargestCoordinate() > 1E-11) {
			Console.WriteLine("versorInverse() test failed (largestCoordinate of VVI = " + (double)Y.LargestCoordinate() + ")");
			return 0; // failure
		}
		
	}
	return 1; // success
}

static int test_versorInverse_dont_mangle_639(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	e3_t V;
	e3_t VI;
	mv gmvV, gmvVI, VVI, VIV, X, Y;
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		V = random_e3_t_dont_mangle_75(genrand() + 0.5);
		
		// avoid 'near'-singular versors
		if (V.LargestCoordinate() > 2.0)
			continue;
		
		// compute inverse
		VI = versorInverse(V);
		
		// convert to GMV
		gmvV = new mv(V);
		gmvVI = new mv(VI);
		
		// compute (inverse(V) V) and (V inverse(V))
		VIV = gp(gmvVI, gmvV);
		VVI = gp(gmvV, gmvVI);
		
		// check that scalar parts are 1
		n1 = VIV.get_scalar();
		n2 = VVI.get_scalar();
		
		if (Math.Abs(n1 - 1.0) > 1E-11) {
			Console.WriteLine("versorInverse() test failed |VI . V - 1|= " + (double)Math.Abs(n1 - 1.0) + ")");
			return 0; // failure
		}
		
		if (Math.Abs(n2 - 1.0) > 1E-11) {
			Console.WriteLine("versorInverse() test failed ( |V . VI - 1| = " + (double)Math.Abs(n2 - 1.0) + ")");
			return 0; // failure
		}
		
		// check that other grade parts are zero:
		X = extractGrade(VIV, GroupBitmap.ALL_GROUPS ^ GroupBitmap.GROUP_0);
		Y = extractGrade(VVI, GroupBitmap.ALL_GROUPS ^ GroupBitmap.GROUP_0);
		
		if (X.LargestCoordinate() > 1E-11) {
			Console.WriteLine("versorInverse() test failed (largestCoordinate of VIV = " + (double)X.LargestCoordinate() + ")");
			return 0; // failure
		}
		
		if (Y.LargestCoordinate() > 1E-11) {
			Console.WriteLine("versorInverse() test failed (largestCoordinate of VVI = " + (double)Y.LargestCoordinate() + ")");
			return 0; // failure
		}
		
	}
	return 1; // success
}

static int test_versorInverse_dont_mangle_640(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 1;
	I3_t V;
	trivector VI;
	mv gmvV, gmvVI, VVI, VIV, X, Y;
	int i;
	double n1, n2;
	
	for (i = 0; i < NB_LOOPS; i++) {
		// get random blade
		V = random_I3_t_dont_mangle_13(genrand() + 0.5);
		
		// avoid 'near'-singular versors
		if (V.LargestCoordinate() > 2.0)
			continue;
		
		// compute inverse
		VI = versorInverse(V);
		
		// convert to GMV
		gmvV = new mv(V);
		gmvVI = new mv(VI);
		
		// compute (inverse(V) V) and (V inverse(V))
		VIV = gp(gmvVI, gmvV);
		VVI = gp(gmvV, gmvVI);
		
		// check that scalar parts are 1
		n1 = VIV.get_scalar();
		n2 = VVI.get_scalar();
		
		if (Math.Abs(n1 - 1.0) > 1E-11) {
			Console.WriteLine("versorInverse() test failed |VI . V - 1|= " + (double)Math.Abs(n1 - 1.0) + ")");
			return 0; // failure
		}
		
		if (Math.Abs(n2 - 1.0) > 1E-11) {
			Console.WriteLine("versorInverse() test failed ( |V . VI - 1| = " + (double)Math.Abs(n2 - 1.0) + ")");
			return 0; // failure
		}
		
		// check that other grade parts are zero:
		X = extractGrade(VIV, GroupBitmap.ALL_GROUPS ^ GroupBitmap.GROUP_0);
		Y = extractGrade(VVI, GroupBitmap.ALL_GROUPS ^ GroupBitmap.GROUP_0);
		
		if (X.LargestCoordinate() > 1E-11) {
			Console.WriteLine("versorInverse() test failed (largestCoordinate of VIV = " + (double)X.LargestCoordinate() + ")");
			return 0; // failure
		}
		
		if (Y.LargestCoordinate() > 1E-11) {
			Console.WriteLine("versorInverse() test failed (largestCoordinate of VVI = " + (double)Y.LargestCoordinate() + ")");
			return 0; // failure
		}
		
	}
	return 1; // success
}

static int test_zero_dont_mangle_641(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 8;
	mv A;
	double s, eps = 0.2;
	int g, i;
	bool z, l, c;
	int basisVectorBitmap = -1;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		g = (int)(genrand() * 3.5);
		A = random_versor_dont_mangle_1_returns_mv(s, g, basisVectorBitmap);
		
		// ask if zero thinks A is zero
		z = zero(A, eps);
		
		// check it yourself
		c = true; // assume it is zero
		for (g = 0; g < 4; g++) {
			if (A.m_c[g] != null) {
				for (int e = 0; e < A.m_c[g].Length; e++) {
					if (Math.Abs(A.m_c[g][e]) > eps) c = false;
				}
			}
		}
		
		// also use mv_largestCoordinate() to verify
		l = !(A.LargestCoordinate() > eps);
		
		if ((c != z) || (l != z)) {
			Console.WriteLine("zero() test failed\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_zero_dont_mangle_642(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	vector A;
	double s, eps = 0.2;
	int i;
	bool l, z;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_vector_dont_mangle_2(s);
		
		// ask if zero thinks A is zero
		z = zero(A, eps);
		
		// use vector_largestCoordinate() to verify
		l = !(A.LargestCoordinate() > eps);
		
		if (l != z) {
			Console.WriteLine("zero() test failed\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_zero_dont_mangle_643(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 4;
	bivector A;
	double s, eps = 0.2;
	int i;
	bool l, z;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_bivector_dont_mangle_4(s);
		
		// ask if zero thinks A is zero
		z = zero(A, eps);
		
		// use bivector_largestCoordinate() to verify
		l = !(A.LargestCoordinate() > eps);
		
		if (l != z) {
			Console.WriteLine("zero() test failed\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_zero_dont_mangle_644(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 2;
	trivector A;
	double s, eps = 0.2;
	int i;
	bool l, z;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_trivector_dont_mangle_7(s);
		
		// ask if zero thinks A is zero
		z = zero(A, eps);
		
		// use trivector_largestCoordinate() to verify
		l = !(A.LargestCoordinate() > eps);
		
		if (l != z) {
			Console.WriteLine("zero() test failed\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_zero_dont_mangle_645(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 5;
	rotor A;
	double s, eps = 0.2;
	int i;
	bool l, z;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_rotor_dont_mangle_8(s);
		
		// ask if zero thinks A is zero
		z = zero(A, eps);
		
		// use rotor_largestCoordinate() to verify
		l = !(A.LargestCoordinate() > eps);
		
		if (l != z) {
			Console.WriteLine("zero() test failed\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_zero_dont_mangle_646(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 2;
	I3_t A;
	double s, eps = 0.2;
	int i;
	bool l, z;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_I3_t_dont_mangle_13(s);
		
		// ask if zero thinks A is zero
		z = zero(A, eps);
		
		// use I3_t_largestCoordinate() to verify
		l = !(A.LargestCoordinate() > eps);
		
		if (l != z) {
			Console.WriteLine("zero() test failed\n");
			return 0; // failure
		}
	}
	return 1; // success
}

static int test_zero_dont_mangle_647(int NB_TESTS_SCALER) 
{
	int NB_LOOPS = 100 + NB_TESTS_SCALER / 2;
	e1_t A;
	double s, eps = 0.2;
	int i;
	bool l, z;
	
	for (i = 0; i < NB_LOOPS; i++) {
		s = genrand();
		A = random_e1_t_dont_mangle_10(s);
		
		// ask if zero thinks A is zero
		z = zero(A, eps);
		
		// use e1_t_largestCoordinate() to verify
		l = !(A.LargestCoordinate() > eps);
		
		if (l != z) {
			Console.WriteLine("zero() test failed\n");
			return 0; // failure
		}
	}
	return 1; // success
}

public static void Main(String[] args)
{
	int retVal = 0;
	// The number of tests will be proportional to this value.
	// This should become a command-line option
	int NB_TESTS_SCALER = 16384;
	
	// seed random number generators with current time
	c3ga.genrand_timeSeed();

	// run all test functions
	if (test_metric_default_mv(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_parse_mv(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_genrand_double(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_add_dont_mangle_382(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_add_dont_mangle_387(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_add_dont_mangle_385(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_add_dont_mangle_388(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_add_dont_mangle_384(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_add_dont_mangle_386(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_add_dont_mangle_383(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_subtract_dont_mangle_389(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_subtract_dont_mangle_390(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_subtract_dont_mangle_391(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_subtract_dont_mangle_392(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_subtract_dont_mangle_393(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyVersor_dont_mangle_394(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyUnitVersor_dont_mangle_395(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyVersorWI_dont_mangle_396(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyVersor_dont_mangle_397(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyUnitVersor_dont_mangle_398(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyVersorWI_dont_mangle_399(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyVersor_dont_mangle_400(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyUnitVersor_dont_mangle_401(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyVersorWI_dont_mangle_402(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyVersor_dont_mangle_403(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyUnitVersor_dont_mangle_404(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyVersorWI_dont_mangle_405(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyVersor_dont_mangle_406(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyUnitVersor_dont_mangle_407(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyVersor_dont_mangle_408(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyUnitVersor_dont_mangle_409(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyVersorWI_dont_mangle_410(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyOM_dont_mangle_411(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyOM_dont_mangle_412(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyOM_dont_mangle_413(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyOM_dont_mangle_414(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyOM_dont_mangle_415(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyOM_dont_mangle_416(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyOM_dont_mangle_417(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyOM_dont_mangle_418(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyOM_dont_mangle_419(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyOM_dont_mangle_420(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_applyOM_dont_mangle_421(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_div_dont_mangle_422(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_div_dont_mangle_423(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_div_dont_mangle_424(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_div_dont_mangle_425(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_div_dont_mangle_426(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_div_dont_mangle_427(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_div_dont_mangle_428(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_dual_dont_mangle_429(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_undual_dont_mangle_430(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_dual_dont_mangle_433(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_undual_dont_mangle_434(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_dual_dont_mangle_435(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_undual_dont_mangle_436(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_dual_dont_mangle_437(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_undual_dont_mangle_438(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_dual_dont_mangle_439(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_undual_dont_mangle_440(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_dual_dont_mangle_441(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_undual_dont_mangle_442(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_dual_dont_mangle_443(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_undual_dont_mangle_444(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_dual_dont_mangle_445(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_undual_dont_mangle_446(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_equals_dont_mangle_447(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_equals_dont_mangle_448(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_equals_dont_mangle_449(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_equals_dont_mangle_450(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_equals_dont_mangle_451(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_equals_dont_mangle_452(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_equals_dont_mangle_453(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_equals_dont_mangle_454(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_equals_dont_mangle_455(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_extractGrade_dont_mangle_456(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_extractGrade2_dont_mangle_457(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_extractGrade0_dont_mangle_458(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_extractGrade2_dont_mangle_459(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_extractGrade0_dont_mangle_460(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_extractGrade1_dont_mangle_461(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_extractGrade2_dont_mangle_462(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_extractGrade3_dont_mangle_463(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_extractGrade0_dont_mangle_464(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_extractGrade1_dont_mangle_465(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_extractGrade2_dont_mangle_466(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_extractGrade3_dont_mangle_467(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_extractGrade0_dont_mangle_468(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_extractGrade1_dont_mangle_469(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_extractGrade2_dont_mangle_470(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_extractGrade3_dont_mangle_471(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gp_dont_mangle_472(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gp_dont_mangle_477(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gp_dont_mangle_475(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gp_dont_mangle_478(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gp_dont_mangle_479(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gp_dont_mangle_474(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gp_dont_mangle_476(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gp_dont_mangle_473(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gp_dont_mangle_480(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gradeBitmap_dont_mangle_481(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gradeBitmap_dont_mangle_482(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gradeBitmap_dont_mangle_483(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gradeBitmap_dont_mangle_484(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gradeBitmap_dont_mangle_485(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gradeBitmap_dont_mangle_486(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gradeBitmap_dont_mangle_487(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gradeBitmap_dont_mangle_488(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hp_dont_mangle_489(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hp_dont_mangle_490(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hp_dont_mangle_491(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hp_dont_mangle_492(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hp_dont_mangle_493(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hp_dont_mangle_494(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hp_dont_mangle_495(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hp_dont_mangle_496(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hp_dont_mangle_497(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hp_dont_mangle_498(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hp_dont_mangle_499(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_ihp_dont_mangle_500(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_ihp_dont_mangle_501(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_ihp_dont_mangle_502(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_ihp_dont_mangle_503(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_ihp_dont_mangle_504(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_ihp_dont_mangle_505(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_ihp_dont_mangle_506(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_ihp_dont_mangle_507(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_ihp_dont_mangle_508(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_ihp_dont_mangle_509(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_increment_dont_mangle_510(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_increment_dont_mangle_511(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_increment_dont_mangle_512(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_decrement_dont_mangle_513(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_decrement_dont_mangle_514(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_decrement_dont_mangle_515(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_sp_dont_mangle_516(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_lc_dont_mangle_517(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_rc_dont_mangle_518(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hip_dont_mangle_519(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_mhip_dont_mangle_520(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_sp_dont_mangle_521(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_lc_dont_mangle_522(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_rc_dont_mangle_523(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hip_dont_mangle_524(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_mhip_dont_mangle_525(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_sp_dont_mangle_526(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_lc_dont_mangle_527(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_rc_dont_mangle_528(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hip_dont_mangle_529(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_mhip_dont_mangle_530(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_sp_dont_mangle_531(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_lc_dont_mangle_532(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_rc_dont_mangle_533(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hip_dont_mangle_534(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_mhip_dont_mangle_535(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_sp_dont_mangle_536(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_lc_dont_mangle_537(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_rc_dont_mangle_538(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hip_dont_mangle_539(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_mhip_dont_mangle_540(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_sp_dont_mangle_541(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_lc_dont_mangle_542(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_rc_dont_mangle_543(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hip_dont_mangle_544(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_mhip_dont_mangle_545(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_sp_dont_mangle_546(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_lc_dont_mangle_547(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_rc_dont_mangle_548(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hip_dont_mangle_549(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_mhip_dont_mangle_550(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_sp_dont_mangle_551(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_lc_dont_mangle_552(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_rc_dont_mangle_553(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hip_dont_mangle_554(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_mhip_dont_mangle_555(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_sp_dont_mangle_556(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_lc_dont_mangle_557(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_rc_dont_mangle_558(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_hip_dont_mangle_559(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_mhip_dont_mangle_560(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_norm_dont_mangle_561(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_norm_dont_mangle_564(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_norm_dont_mangle_566(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_norm_dont_mangle_567(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_norm_dont_mangle_563(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_norm_dont_mangle_565(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_norm_dont_mangle_562(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_norm_dont_mangle_568(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_norm2_dont_mangle_569(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_norm2_dont_mangle_570(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_norm2_dont_mangle_571(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_norm2_dont_mangle_572(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_norm2_dont_mangle_573(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_norm2_dont_mangle_574(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_norm2_dont_mangle_575(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_norm2_dont_mangle_576(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_op_dont_mangle_577_1(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_op_dont_mangle_577_2(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_op_dont_mangle_578(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_op_dont_mangle_579(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_op_dont_mangle_580(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_op_dont_mangle_581(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_op_dont_mangle_582(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_op_dont_mangle_583(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_op_dont_mangle_584(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_exp_dont_mangle_586(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_exp_dont_mangle_587(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_cosh_dont_mangle_589(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_cosh_dont_mangle_588(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_sinh_dont_mangle_590(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_sinh_dont_mangle_591(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_cos_dont_mangle_593(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_sin_dont_mangle_595(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_negate_dont_mangle_596(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_cliffordConjugate_dont_mangle_597(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gradeInvolution_dont_mangle_598(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_reverse_dont_mangle_599(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_negate_dont_mangle_600(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_cliffordConjugate_dont_mangle_601(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gradeInvolution_dont_mangle_602(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_reverse_dont_mangle_603(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_negate_dont_mangle_604(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_cliffordConjugate_dont_mangle_605(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gradeInvolution_dont_mangle_606(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_reverse_dont_mangle_607(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_negate_dont_mangle_608(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_cliffordConjugate_dont_mangle_609(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gradeInvolution_dont_mangle_610(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_reverse_dont_mangle_611(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_negate_dont_mangle_612(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_cliffordConjugate_dont_mangle_613(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gradeInvolution_dont_mangle_614(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_reverse_dont_mangle_615(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_negate_dont_mangle_616(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_cliffordConjugate_dont_mangle_617(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_gradeInvolution_dont_mangle_618(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_reverse_dont_mangle_619(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_unit_dont_mangle_624(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_unit_dont_mangle_626(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_unit_dont_mangle_628(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_unit_dont_mangle_627(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_unit_dont_mangle_625(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_unit_dont_mangle_629(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_unit_dont_mangle_630(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_unit_dont_mangle_631(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_unit_dont_mangle_632(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_versorInverse_dont_mangle_633(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_versorInverse_dont_mangle_634(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_versorInverse_dont_mangle_635(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_versorInverse_dont_mangle_636(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_versorInverse_dont_mangle_637(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_versorInverse_dont_mangle_638(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_versorInverse_dont_mangle_639(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_versorInverse_dont_mangle_640(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_zero_dont_mangle_641(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_zero_dont_mangle_642(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_zero_dont_mangle_643(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_zero_dont_mangle_644(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_zero_dont_mangle_645(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_zero_dont_mangle_646(NB_TESTS_SCALER) == 0) retVal = -1;
	if (test_zero_dont_mangle_647(NB_TESTS_SCALER) == 0) retVal = -1;

	if (retVal != 0) Console.WriteLine("Test failed.\n");
	else Console.WriteLine("Done.\n");	

	Environment.Exit(retVal);
}
} // end of class TestSuite
} // end of namespace c3ga_ns
