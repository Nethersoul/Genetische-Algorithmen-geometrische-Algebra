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

[FlagsAttribute]
public enum GroupBitmap : int
{
	GROUP_0  = 1, //1
	GROUP_1  = 2, //e1, e2, e3
	GROUP_2  = 4, //e1^e2, e1^e3, e2^e3
	GROUP_3  = 8, //e1^e2^e3

	ALL_GROUPS = 15, // all groups

	GRADE_0 = 1, 
	GRADE_1 = 2, 
	GRADE_2 = 4, 
	GRADE_3 = 8, 

	ALL_GRADES = 15 // all grades
}

/// <summary>
/// These constants define a unique number for each specialized multivector type.
/// They are used to report usage of non-optimized functions.
/// </summary>
public enum SmvType {
	C3GA_NONE = -1,
	C3GA_MV = 0,
	C3GA_DOUBLE = 1,
	C3GA_E1_T = 2,
	C3GA_E2_T = 3,
	C3GA_E3_T = 4,
	C3GA_I3_T = 5,
	C3GA_VECTOR = 6,
	C3GA_BIVECTOR = 7,
	C3GA_TRIVECTOR = 8,
	C3GA_ROTOR = 9,
	C3GA_ODDVERSOR = 10,
	C3GA_INVALID
}

public class c3ga 
{ 
	public static readonly e1_t e1 = new e1_t();
	public static readonly e2_t e2 = new e2_t();
	public static readonly e3_t e3 = new e3_t();
	public static readonly I3_t I3 = new I3_t();
	public static readonly vector vectorE1 = new vector(vector.coord_e1_e2_e3, 1.0, 0.0, 0.0);
	public static readonly vector vectorE2 = new vector(vector.coord_e1_e2_e3, 0.0, 2.0, 0.0);
	public static readonly vector vectorE3 = new vector(vector.coord_e1_e2_e3, 0.0, 0.0, 3.0);
	public static readonly bivector someBivectorConstant = new bivector(bivector.coord_e1e2_e2e3_e3e1, 1.0, 2.0, 3.0);
	public static readonly rotor rotor90 = new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1, 0.707111, 0.707111, 0.0, 0.0);


	static protected internal string[] typenames = 
		new string[] {
			"mv",
			"double",
			"e1_t",
			"e2_t",
			"e3_t",
			"I3_t",
			"vector",
			"bivector",
			"trivector",
			"rotor",
			"oddVersor"
		};

	/// <summary>The dimension of the space
	/// </summary>
	public const int SpaceDim = 3;
	/// <summary>Number of groups/grades of coordinates in a multivector
	/// </summary>
	public const int NbGroups = 4;
	/// <summary>Is the metric of the space Euclidean? (false or true)
	/// </summary>
	public const bool MetricEuclidean = true;
	/// <summary>Names of the basis vectors.
	/// </summary>
	public static readonly string[] BasisVectorNames = new string[] {
		"e1", "e2", "e3"
	};
	/// <summary>The constants for the grades, in an array.
	/// </summary>
	public static readonly GroupBitmap[] Grades = {GroupBitmap.GRADE_0, GroupBitmap.GRADE_1, GroupBitmap.GRADE_2, GroupBitmap.GRADE_3, 0, 0, 0, 0};
	/// <summary>The constants for the groups, in an array.
	/// </summary>
	public static readonly GroupBitmap[] Groups = {GroupBitmap.GROUP_0, GroupBitmap.GROUP_1, GroupBitmap.GROUP_2, GroupBitmap.GROUP_3};
	/// <summary>This array can be used to lookup the number of coordinates for a group part of a general multivector.
	/// </summary>
	public static readonly int[] GroupSize = { 1, 3, 3, 1 };
	/// <summary>This array can be used to lookup the number of coordinates based on a group usage bitmap.
	/// </summary>
	public static readonly int[] MvSize = new int[] {
		0, 1, 3, 4, 3, 4, 6, 7, 1, 2, 4, 5, 4, 5, 7, 8	};
	/// <summary>This array of integers contains the 'sign' (even/odd permutation of canonical order) of basis elements in the general multivector.
	/// Use it to answer 'what is the permutation of the coordinate at index [x]'?
	/// </summary>
	public static readonly double[] BasisElementSignByIndex = new double[]
		{1, 1, 1, 1, 1, 1, 1, 1};
	/// <summary>This array of integers contains the 'sign' (even/odd permutation of canonical order) of basis elements in the general multivector.
	/// Use it to answer 'what is the permutation of the coordinate of bitmap [x]'?
	/// </summary>
	public static readonly double[] BasisElementSignByBitmap = new double[]
		{1, 1, 1, 1, 1, 1, 1, 1};
	/// <summary>This array of integers contains the order of basis elements in the general multivector.
	/// Use it to answer: 'at what index do I find basis element [x] (x = basis vector bitmap)?'
	/// </summary>
	public static readonly int[] BasisElementIndexByBitmap = new int[]
		{0, 1, 2, 4, 3, 5, 6, 7};
	/// <summary>This array of integers contains the indices of basis elements in the general multivector.
	/// Use it to answer: 'what basis element do I find at index [x]'?
	/// </summary>
	public static readonly int[] BasisElementBitmapByIndex = new int[]
		{0, 1, 2, 4, 3, 5, 6, 7};
	/// <summary>This array of grade of each basis elements in the general multivector.
	/// Use it to answer: 'what is the grade of basis element bitmap [x]'?
	/// </summary>
	public static readonly int[] BasisElementGradeByBitmap = new int[]
		{0, 1, 1, 2, 1, 2, 2, 3};
	/// <summary>This array of group of each basis elements in the general multivector.
	/// Use it to answer: 'what is the group of basis element bitmap [x]'?
	/// </summary>
	public static readonly int[] BasisElementGroupByBitmap = new int[]
		{0, 1, 1, 2, 1, 2, 2, 3};
	/// <summary>This array of integers contains the order of basis elements in the general multivector.
	/// Use it to answer: 'what basis vectors are in the basis element at position [x]?
	/// </summary>
	public static readonly int[][] BasisElements = new int[][] {
		new int[] {-1},
		new int[] {0, -1},
		new int[] {1, -1},
		new int[] {2, -1},
		new int[] {0, 1, -1},
		new int[] {0, 2, -1},
		new int[] {1, 2, -1},
		new int[] {0, 1, 2, -1}
	};

    // Simple code for multi-threading-safe random number generation (the Random class is not thread-safe).
    // The worst that can happen is that multiple threads get the same random number.
    protected internal static System.Random s_randomGenerator = new System.Random();
    protected internal const int NB_RANDOM_DOUBLES = 32; // must be power of two
    protected internal static double[] s_randomDoubles = new double[NB_RANDOM_DOUBLES];
    protected internal static int s_randomDoublesIdx = NB_RANDOM_DOUBLES;
    protected internal static double NextRandomDouble() {
        if (s_randomDoublesIdx >= s_randomDoubles.Length) {
            lock(s_randomDoubles) {
                for (int i = 0; i < s_randomDoubles.Length; i++) {
                    s_randomDoubles[i] = s_randomGenerator.NextDouble();
                }
                s_randomDoublesIdx = 0;
            }
        }
        int idx = s_randomDoublesIdx & (NB_RANDOM_DOUBLES-1);
        s_randomDoublesIdx++;
        return s_randomDoubles[idx];
	}

    /// <summary>Sets 1 doubles to zero.</summary>
	protected internal static void Zero_1(double[] dst) {
		dst[0]=0.0;
	}
	/// <summary>Copies 1 doubles from 'src' to 'dst'. </summary>
	protected internal static void Copy_1(double[] dst, double[] src) {
			dst[0] = src[0];
	}
    /// <summary>Sets 2 doubles to zero.</summary>
	protected internal static void Zero_2(double[] dst) {
		dst[0]=dst[1]=0.0;
	}
	/// <summary>Copies 2 doubles from 'src' to 'dst'. </summary>
	protected internal static void Copy_2(double[] dst, double[] src) {
			dst[0] = src[0];
			dst[1] = src[1];
	}
    /// <summary>Sets 3 doubles to zero.</summary>
	protected internal static void Zero_3(double[] dst) {
		dst[0]=dst[1]=dst[2]=0.0;
	}
	/// <summary>Copies 3 doubles from 'src' to 'dst'. </summary>
	protected internal static void Copy_3(double[] dst, double[] src) {
			dst[0] = src[0];
			dst[1] = src[1];
			dst[2] = src[2];
	}
    /// <summary>Sets 4 doubles to zero.</summary>
	protected internal static void Zero_4(double[] dst) {
		dst[0]=dst[1]=dst[2]=dst[3]=0.0;
	}
	/// <summary>Copies 4 doubles from 'src' to 'dst'. </summary>
	protected internal static void Copy_4(double[] dst, double[] src) {
			dst[0] = src[0];
			dst[1] = src[1];
			dst[2] = src[2];
			dst[3] = src[3];
	}
    /// <summary>Sets 5 doubles to zero.</summary>
	protected internal static void Zero_5(double[] dst) {
		dst[0]=dst[1]=dst[2]=dst[3]=dst[4]=0.0;
	}
	/// <summary>Copies 5 doubles from 'src' to 'dst'. </summary>
	protected internal static void Copy_5(double[] dst, double[] src) {
			dst[0] = src[0];
			dst[1] = src[1];
			dst[2] = src[2];
			dst[3] = src[3];
			dst[4] = src[4];
	}
    /// <summary>Sets 6 doubles to zero.</summary>
	protected internal static void Zero_6(double[] dst) {
		dst[0]=dst[1]=dst[2]=dst[3]=dst[4]=dst[5]=0.0;
	}
	/// <summary>Copies 6 doubles from 'src' to 'dst'. </summary>
	protected internal static void Copy_6(double[] dst, double[] src) {
			dst[0] = src[0];
			dst[1] = src[1];
			dst[2] = src[2];
			dst[3] = src[3];
			dst[4] = src[4];
			dst[5] = src[5];
	}
    /// <summary>Sets 7 doubles to zero.</summary>
	protected internal static void Zero_7(double[] dst) {
		dst[0]=dst[1]=dst[2]=dst[3]=dst[4]=dst[5]=dst[6]=0.0;
	}
	/// <summary>Copies 7 doubles from 'src' to 'dst'. </summary>
	protected internal static void Copy_7(double[] dst, double[] src) {
			dst[0] = src[0];
			dst[1] = src[1];
			dst[2] = src[2];
			dst[3] = src[3];
			dst[4] = src[4];
			dst[5] = src[5];
			dst[6] = src[6];
	}
    /// <summary>Sets 8 doubles to zero.</summary>
	protected internal static void Zero_8(double[] dst) {
		dst[0]=dst[1]=dst[2]=dst[3]=dst[4]=dst[5]=dst[6]=dst[7]=0.0;
	}
	/// <summary>Copies 8 doubles from 'src' to 'dst'. </summary>
	protected internal static void Copy_8(double[] dst, double[] src) {
			dst[0] = src[0];
			dst[1] = src[1];
			dst[2] = src[2];
			dst[3] = src[3];
			dst[4] = src[4];
			dst[5] = src[5];
			dst[6] = src[6];
			dst[7] = src[7];
	}
    /// <summary>Sets 9 doubles to zero.</summary>
	protected internal static void Zero_9(double[] dst) {
		dst[0]=dst[1]=dst[2]=dst[3]=dst[4]=dst[5]=dst[6]=dst[7]=dst[8]=0.0;
	}
	/// <summary>Copies 9 doubles from 'src' to 'dst'. </summary>
	protected internal static void Copy_9(double[] dst, double[] src) {
			dst[0] = src[0];
			dst[1] = src[1];
			dst[2] = src[2];
			dst[3] = src[3];
			dst[4] = src[4];
			dst[5] = src[5];
			dst[6] = src[6];
			dst[7] = src[7];
			dst[8] = src[8];
	}
    /// <summary>Sets 10 doubles to zero.</summary>
	protected internal static void Zero_10(double[] dst) {
		dst[0]=dst[1]=dst[2]=dst[3]=dst[4]=dst[5]=dst[6]=dst[7]=dst[8]=dst[9]=0.0;
	}
	/// <summary>Copies 10 doubles from 'src' to 'dst'. </summary>
	protected internal static void Copy_10(double[] dst, double[] src) {
			dst[0] = src[0];
			dst[1] = src[1];
			dst[2] = src[2];
			dst[3] = src[3];
			dst[4] = src[4];
			dst[5] = src[5];
			dst[6] = src[6];
			dst[7] = src[7];
			dst[8] = src[8];
			dst[9] = src[9];
	}
    /// <summary>Sets 11 doubles to zero.</summary>
	protected internal static void Zero_11(double[] dst) {
		dst[0]=dst[1]=dst[2]=dst[3]=dst[4]=dst[5]=dst[6]=dst[7]=dst[8]=dst[9]=dst[10]=0.0;
	}
	/// <summary>Copies 11 doubles from 'src' to 'dst'. </summary>
	protected internal static void Copy_11(double[] dst, double[] src) {
			dst[0] = src[0];
			dst[1] = src[1];
			dst[2] = src[2];
			dst[3] = src[3];
			dst[4] = src[4];
			dst[5] = src[5];
			dst[6] = src[6];
			dst[7] = src[7];
			dst[8] = src[8];
			dst[9] = src[9];
			dst[10] = src[10];
	}
    /// <summary>Sets 12 doubles to zero.</summary>
	protected internal static void Zero_12(double[] dst) {
		dst[0]=dst[1]=dst[2]=dst[3]=dst[4]=dst[5]=dst[6]=dst[7]=dst[8]=dst[9]=dst[10]=dst[11]=0.0;
	}
	/// <summary>Copies 12 doubles from 'src' to 'dst'. </summary>
	protected internal static void Copy_12(double[] dst, double[] src) {
			dst[0] = src[0];
			dst[1] = src[1];
			dst[2] = src[2];
			dst[3] = src[3];
			dst[4] = src[4];
			dst[5] = src[5];
			dst[6] = src[6];
			dst[7] = src[7];
			dst[8] = src[8];
			dst[9] = src[9];
			dst[10] = src[10];
			dst[11] = src[11];
	}
    /// <summary>Sets 13 doubles to zero.</summary>
	protected internal static void Zero_13(double[] dst) {
		dst[0]=dst[1]=dst[2]=dst[3]=dst[4]=dst[5]=dst[6]=dst[7]=dst[8]=dst[9]=dst[10]=dst[11]=dst[12]=0.0;
	}
	/// <summary>Copies 13 doubles from 'src' to 'dst'. </summary>
	protected internal static void Copy_13(double[] dst, double[] src) {
			dst[0] = src[0];
			dst[1] = src[1];
			dst[2] = src[2];
			dst[3] = src[3];
			dst[4] = src[4];
			dst[5] = src[5];
			dst[6] = src[6];
			dst[7] = src[7];
			dst[8] = src[8];
			dst[9] = src[9];
			dst[10] = src[10];
			dst[11] = src[11];
			dst[12] = src[12];
	}
    /// <summary>Sets 14 doubles to zero.</summary>
	protected internal static void Zero_14(double[] dst) {
		dst[0]=dst[1]=dst[2]=dst[3]=dst[4]=dst[5]=dst[6]=dst[7]=dst[8]=dst[9]=dst[10]=dst[11]=dst[12]=dst[13]=0.0;
	}
	/// <summary>Copies 14 doubles from 'src' to 'dst'. </summary>
	protected internal static void Copy_14(double[] dst, double[] src) {
			dst[0] = src[0];
			dst[1] = src[1];
			dst[2] = src[2];
			dst[3] = src[3];
			dst[4] = src[4];
			dst[5] = src[5];
			dst[6] = src[6];
			dst[7] = src[7];
			dst[8] = src[8];
			dst[9] = src[9];
			dst[10] = src[10];
			dst[11] = src[11];
			dst[12] = src[12];
			dst[13] = src[13];
	}
    /// <summary>Sets 15 doubles to zero.</summary>
	protected internal static void Zero_15(double[] dst) {
		dst[0]=dst[1]=dst[2]=dst[3]=dst[4]=dst[5]=dst[6]=dst[7]=dst[8]=dst[9]=dst[10]=dst[11]=dst[12]=dst[13]=dst[14]=0.0;
	}
	/// <summary>Copies 15 doubles from 'src' to 'dst'. </summary>
	protected internal static void Copy_15(double[] dst, double[] src) {
			dst[0] = src[0];
			dst[1] = src[1];
			dst[2] = src[2];
			dst[3] = src[3];
			dst[4] = src[4];
			dst[5] = src[5];
			dst[6] = src[6];
			dst[7] = src[7];
			dst[8] = src[8];
			dst[9] = src[9];
			dst[10] = src[10];
			dst[11] = src[11];
			dst[12] = src[12];
			dst[13] = src[13];
			dst[14] = src[14];
	}
    /// <summary>Sets 16 doubles to zero.</summary>
	protected internal static void Zero_16(double[] dst) {
		dst[0]=dst[1]=dst[2]=dst[3]=dst[4]=dst[5]=dst[6]=dst[7]=dst[8]=dst[9]=dst[10]=dst[11]=dst[12]=dst[13]=dst[14]=dst[15]=0.0;
	}
	/// <summary>Copies 16 doubles from 'src' to 'dst'. </summary>
	protected internal static void Copy_16(double[] dst, double[] src) {
			dst[0] = src[0];
			dst[1] = src[1];
			dst[2] = src[2];
			dst[3] = src[3];
			dst[4] = src[4];
			dst[5] = src[5];
			dst[6] = src[6];
			dst[7] = src[7];
			dst[8] = src[8];
			dst[9] = src[9];
			dst[10] = src[10];
			dst[11] = src[11];
			dst[12] = src[12];
			dst[13] = src[13];
			dst[14] = src[14];
			dst[15] = src[15];
	}
	/// <summary>Sets N doubles to zero.</summary>
	protected internal static void Zero_N(double[] dst, int N) {
		for (int i = 0; i < N; i++)
			dst[i] = 0.0;
	}
	/// <summary>Copies N doubles from 'src' to 'dst'. </summary>
	protected internal static void Copy_N(double[] dst, double[] src, int N) {
		for (int i = 0; i < N; i++)
			dst[i] = src[i];
	}


	private const string DEFAULT_FP = "F2";
	private const string DEFAULT_START = "";
	private const string DEFAULT_END = "";
	private const string DEFAULT_MUL = "*";
	private const string DEFAULT_WEDGE = "^";
	private const string DEFAULT_PLUS = " + ";
	private const string DEFAULT_MINUS = " - ";

	// These strings determine how the output of string() is formatted.
	// You can alter them at runtime using setStringFormat(). 
	protected static string string_fp = DEFAULT_FP;
	protected static string string_start = DEFAULT_START;
	protected static string string_end = DEFAULT_END;
	protected static string string_mul = DEFAULT_MUL;
	protected static string string_wedge = DEFAULT_WEDGE;
	protected static string string_plus = DEFAULT_PLUS;
	protected static string string_minus = DEFAULT_MINUS;
	
	public const string STRING_FP = "fp";
	public const string STRING_START = "start";
	public const string STRING_END = "end";
	public const string STRING_MUL = "mul";
	public const string STRING_WEDGE = "wedge";
	public const string STRING_PLUS = "plus";
	public const string STRING_MINUS= "minus";

	/// <summary>
	/// Sets the formatting of toString().
	/// </summary>
	/// 
	/// <param name="what">What formatter to set. Valid values: string_FP, string_START, string_END, string_MUL, string_WEDGE, string_PLUS, string_MINUS.</param>
	/// <param name="format">The value for 'what'. Use 'null' to set the default value.</param>
	public static void SetStringFormat(string what, string format) {
		if (what.Equals(STRING_FP)) 
			string_fp = (format != null) ? format : DEFAULT_FP;
		else if (what.Equals(STRING_START)) 
			string_start = (format != null) ? format : DEFAULT_START;
		else if (what.Equals(STRING_END)) 
			string_end = (format != null) ? format : DEFAULT_END;
		else if (what.Equals(STRING_MUL)) 
			string_mul = (format != null) ? format : DEFAULT_MUL;
		else if (what.Equals(STRING_WEDGE)) 
			string_wedge = (format != null) ? format : DEFAULT_WEDGE;
		else if (what.Equals(STRING_PLUS)) 
			string_plus = (format != null) ? format : DEFAULT_PLUS;
		else if (what.Equals(STRING_MINUS)) 
			string_minus = (format != null) ? format : DEFAULT_MINUS;
		else throw new ArgumentException("invalid argument to setStringFormat(): " + what);
	}
	
    /// <summary>Converts a multivector to a string using default float format.</summary>
	public static string String(mv_if value) {
		return String(value, null);
	}
	
	/// <summary>
	/// Converts a multivector to a string according to a float format like "F".
	/// </summary>
	/// <param name="fp">floating point format. Use 'null' for the default format (see setStringFormat()).</param>
	public static string String(mv_if value, string fp) {
		mv obj = value.to_mv();
		System.Text.StringBuilder result = new System.Text.StringBuilder();
		int ia = 0; // global index into coordinates (runs from 0 to 7)
		int cnt = 0; // how many coordinates printed so far

		// set up the floating point precision
		if (fp == null) fp = string_fp;

		// start the string
		result.Append(string_start);

		// print all coordinates
		for (int g = 0; g < 4; g++) {
			double[] Cg = obj.m_c[g];
			if (Cg != null) {
				for (int b = 0; b < GroupSize[g]; b++) {
					double coord = (double)BasisElementSignByIndex[ia] * Cg[b];
					
					// goal: print [+|-]obj.m_c[k][* basisVector1 ^ ... ^ basisVectorN]
					
					string tmpFloatStr = Math.Abs(coord).ToString(fp);

					if (Double.Parse(tmpFloatStr) != 0.0) {
						// print [+|-]
						result.Append((coord >= 0.0) 
							? ((cnt>0) ? string_plus : "")
							: string_minus);
						// print obj.m_c[k]
						result.Append(tmpFloatStr);

						if (g != 0) { // if not grade 0, print [* basisVector1 ^ ... ^ basisVectorN]
							result.Append(string_mul);

							// print all basis vectors
							int bei = 0;
							while (BasisElements[ia][bei] >= 0) {
								if (bei > 0)
									result.Append(string_wedge);
								result.Append(BasisVectorNames[BasisElements[ia][bei]]);
								bei++;
							}
						}

						cnt++;
					}
					ia++;
				}
			}
			else ia += GroupSize[g];
		}

		// if no coordinates printed: 0
		if (cnt == 0) result.Append("0");

		// end the string
		result.Append(string_end);

		return result.ToString();
	}


	/// <summary>
    /// Simple way to call parser (regardless of whether it is a builtin or ANTLR parser).
    /// 
    /// Throws a ParseException on failure.
    /// 
    /// When an ANTLR based parser throws an exception, 
    /// all details (like line number and cause) are lost. 
    /// If these details are required, call the ANTLR parser directly.
    /// 
    /// </summary>
    /// <param name="str">The multivector string to be parsed (can be output of mv.ToString())</param>
    /// <returns>Multivector value represented by 'str'.</returns>
    public static mv Parse(string str)
    {
        return Parser.Parse(str, "string");
    }
	/// <summary>Computes the partial geometric product of two multivectors (group 0  x  group 0 -> group 0)
	/// </summary>
	protected internal static void gp_default_0_0_0(double[] A, double[] B, double[] C) {
		C[0] += A[0]*B[0];
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 0  x  group 1 -> group 1)
	/// </summary>
	protected internal static void gp_default_0_1_1(double[] A, double[] B, double[] C) {
		C[0] += A[0]*B[0];
		C[1] += A[0]*B[1];
		C[2] += A[0]*B[2];
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 0  x  group 2 -> group 2)
	/// </summary>
	protected internal static void gp_default_0_2_2(double[] A, double[] B, double[] C) {
		gp_default_0_1_1(A, B, C);
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 0  x  group 3 -> group 3)
	/// </summary>
	protected internal static void gp_default_0_3_3(double[] A, double[] B, double[] C) {
		gp_default_0_0_0(A, B, C);
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 1  x  group 0 -> group 1)
	/// </summary>
	protected internal static void gp_default_1_0_1(double[] A, double[] B, double[] C) {
		C[0] += A[0]*B[0];
		C[1] += A[1]*B[0];
		C[2] += A[2]*B[0];
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 1  x  group 1 -> group 0)
	/// </summary>
	protected internal static void gp_default_1_1_0(double[] A, double[] B, double[] C) {
		C[0] += (A[0]*B[0]+A[1]*B[1]+A[2]*B[2]);
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 1  x  group 1 -> group 2)
	/// </summary>
	protected internal static void gp_default_1_1_2(double[] A, double[] B, double[] C) {
		C[0] += (A[0]*B[1]-A[1]*B[0]);
		C[1] += (A[0]*B[2]-A[2]*B[0]);
		C[2] += (A[1]*B[2]-A[2]*B[1]);
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 1  x  group 2 -> group 1)
	/// </summary>
	protected internal static void gp_default_1_2_1(double[] A, double[] B, double[] C) {
		C[0] += (-A[1]*B[0]-A[2]*B[1]);
		C[1] += (A[0]*B[0]-A[2]*B[2]);
		C[2] += (A[0]*B[1]+A[1]*B[2]);
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 1  x  group 2 -> group 3)
	/// </summary>
	protected internal static void gp_default_1_2_3(double[] A, double[] B, double[] C) {
		C[0] += (A[0]*B[2]-A[1]*B[1]+A[2]*B[0]);
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 1  x  group 3 -> group 2)
	/// </summary>
	protected internal static void gp_default_1_3_2(double[] A, double[] B, double[] C) {
		C[0] += A[2]*B[0];
		C[1] += -A[1]*B[0];
		C[2] += A[0]*B[0];
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 2  x  group 0 -> group 2)
	/// </summary>
	protected internal static void gp_default_2_0_2(double[] A, double[] B, double[] C) {
		gp_default_1_0_1(A, B, C);
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 2  x  group 1 -> group 1)
	/// </summary>
	protected internal static void gp_default_2_1_1(double[] A, double[] B, double[] C) {
		C[0] += (A[0]*B[1]+A[1]*B[2]);
		C[1] += (-A[0]*B[0]+A[2]*B[2]);
		C[2] += (-A[1]*B[0]-A[2]*B[1]);
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 2  x  group 1 -> group 3)
	/// </summary>
	protected internal static void gp_default_2_1_3(double[] A, double[] B, double[] C) {
		gp_default_1_2_3(A, B, C);
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 2  x  group 2 -> group 0)
	/// </summary>
	protected internal static void gp_default_2_2_0(double[] A, double[] B, double[] C) {
		C[0] += (-A[0]*B[0]-A[1]*B[1]-A[2]*B[2]);
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 2  x  group 2 -> group 2)
	/// </summary>
	protected internal static void gp_default_2_2_2(double[] A, double[] B, double[] C) {
		C[0] += (-A[1]*B[2]+A[2]*B[1]);
		C[1] += (A[0]*B[2]-A[2]*B[0]);
		C[2] += (-A[0]*B[1]+A[1]*B[0]);
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 2  x  group 3 -> group 1)
	/// </summary>
	protected internal static void gp_default_2_3_1(double[] A, double[] B, double[] C) {
		C[0] += -A[2]*B[0];
		C[1] += A[1]*B[0];
		C[2] += -A[0]*B[0];
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 3  x  group 0 -> group 3)
	/// </summary>
	protected internal static void gp_default_3_0_3(double[] A, double[] B, double[] C) {
		gp_default_0_0_0(A, B, C);
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 3  x  group 1 -> group 2)
	/// </summary>
	protected internal static void gp_default_3_1_2(double[] A, double[] B, double[] C) {
		C[0] += A[0]*B[2];
		C[1] += -A[0]*B[1];
		C[2] += A[0]*B[0];
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 3  x  group 2 -> group 1)
	/// </summary>
	protected internal static void gp_default_3_2_1(double[] A, double[] B, double[] C) {
		C[0] += -A[0]*B[2];
		C[1] += A[0]*B[1];
		C[2] += -A[0]*B[0];
	}
	/// <summary>Computes the partial geometric product of two multivectors (group 3  x  group 3 -> group 0)
	/// </summary>
	protected internal static void gp_default_3_3_0(double[] A, double[] B, double[] C) {
		C[0] += -A[0]*B[0];
	}
	/// <summary>copies coordinates of group 0
	/// </summary>
	protected internal static void copyGroup_0(double[] A, double[] C) {
		C[0] = A[0];
	}
	/// <summary>copies and multiplies (by s) coordinates of group 0
	/// </summary>
	protected internal static void copyMul_0(double[] A, double[] C, double s) {
		C[0] = A[0]*s;
	}
	/// <summary>copies and divides (by s) coordinates of group 0
	/// </summary>
	protected internal static void copyDiv_0(double[] A, double[] C, double s) {
		C[0] = A[0]/s;
	}
	/// <summary>adds coordinates of group 0 from variable A to C
	/// </summary>
	protected internal static void add_0(double[] A, double[] C) {
		C[0] += A[0];
	}
	/// <summary>subtracts coordinates of group 0 in variable A from C
	/// </summary>
	protected internal static void sub_0(double[] A, double[] C) {
		C[0] -= A[0];
	}
	/// <summary>negate coordinates of group 0 of variable A
	/// </summary>
	protected internal static void neg_0(double[] A, double[] C) {
		C[0] = -A[0];
	}
	/// <summary>adds coordinates of group 0 of variables A and B
	/// </summary>
	protected internal static void add2_0_0(double[] A, double[] B, double[] C) {
		C[0] = (A[0]+B[0]);
	}
	/// <summary>subtracts coordinates of group 0 of variables A from B
	/// </summary>
	protected internal static void sub2_0_0(double[] A, double[] B, double[] C) {
		C[0] = (A[0]-B[0]);
	}
	/// <summary>performs coordinate-wise multiplication of coordinates of group 0 of variables A and B
	/// </summary>
	protected internal static void hp_0_0(double[] A, double[] B, double[] C) {
		C[0] = A[0]*B[0];
	}
	/// <summary>performs coordinate-wise division of coordinates of group 0 of variables A and B
	/// (no checks for divide by zero are made)
	/// </summary>
	protected internal static void ihp_0_0(double[] A, double[] B, double[] C) {
		C[0] = A[0]/((B[0]));
	}
	/// <summary>check for equality up to eps of coordinates of group 0 of variables A and B
	/// </summary>
	protected internal static bool equals_0_0(double[] A, double[] B, double eps) {
		if (((A[0] - B[0]) < -eps) || ((A[0] - B[0]) > eps)) return false;
	return true;
	}
	/// <summary>checks if coordinates of group 0 of variable A are zero up to eps
	/// </summary>
	protected internal static bool zeroGroup_0(double[] A, double eps) {
		if ((A[0] < -eps) || (A[0] > eps)) return false;
		return true;
	}
	/// <summary>copies coordinates of group 1
	/// </summary>
	protected internal static void copyGroup_1(double[] A, double[] C) {
		C[0] = A[0];
		C[1] = A[1];
		C[2] = A[2];
	}
	/// <summary>copies and multiplies (by s) coordinates of group 1
	/// </summary>
	protected internal static void copyMul_1(double[] A, double[] C, double s) {
		C[0] = A[0]*s;
		C[1] = A[1]*s;
		C[2] = A[2]*s;
	}
	/// <summary>copies and divides (by s) coordinates of group 1
	/// </summary>
	protected internal static void copyDiv_1(double[] A, double[] C, double s) {
		C[0] = A[0]/s;
		C[1] = A[1]/s;
		C[2] = A[2]/s;
	}
	/// <summary>adds coordinates of group 1 from variable A to C
	/// </summary>
	protected internal static void add_1(double[] A, double[] C) {
		C[0] += A[0];
		C[1] += A[1];
		C[2] += A[2];
	}
	/// <summary>subtracts coordinates of group 1 in variable A from C
	/// </summary>
	protected internal static void sub_1(double[] A, double[] C) {
		C[0] -= A[0];
		C[1] -= A[1];
		C[2] -= A[2];
	}
	/// <summary>negate coordinates of group 1 of variable A
	/// </summary>
	protected internal static void neg_1(double[] A, double[] C) {
		C[0] = -A[0];
		C[1] = -A[1];
		C[2] = -A[2];
	}
	/// <summary>adds coordinates of group 1 of variables A and B
	/// </summary>
	protected internal static void add2_1_1(double[] A, double[] B, double[] C) {
		C[0] = (A[0]+B[0]);
		C[1] = (A[1]+B[1]);
		C[2] = (A[2]+B[2]);
	}
	/// <summary>subtracts coordinates of group 1 of variables A from B
	/// </summary>
	protected internal static void sub2_1_1(double[] A, double[] B, double[] C) {
		C[0] = (A[0]-B[0]);
		C[1] = (A[1]-B[1]);
		C[2] = (A[2]-B[2]);
	}
	/// <summary>performs coordinate-wise multiplication of coordinates of group 1 of variables A and B
	/// </summary>
	protected internal static void hp_1_1(double[] A, double[] B, double[] C) {
		C[0] = A[0]*B[0];
		C[1] = A[1]*B[1];
		C[2] = A[2]*B[2];
	}
	/// <summary>performs coordinate-wise division of coordinates of group 1 of variables A and B
	/// (no checks for divide by zero are made)
	/// </summary>
	protected internal static void ihp_1_1(double[] A, double[] B, double[] C) {
		C[0] = A[0]/((B[0]));
		C[1] = A[1]/((B[1]));
		C[2] = A[2]/((B[2]));
	}
	/// <summary>check for equality up to eps of coordinates of group 1 of variables A and B
	/// </summary>
	protected internal static bool equals_1_1(double[] A, double[] B, double eps) {
		if (((A[0] - B[0]) < -eps) || ((A[0] - B[0]) > eps)) return false;
		if (((A[1] - B[1]) < -eps) || ((A[1] - B[1]) > eps)) return false;
		if (((A[2] - B[2]) < -eps) || ((A[2] - B[2]) > eps)) return false;
	return true;
	}
	/// <summary>checks if coordinates of group 1 of variable A are zero up to eps
	/// </summary>
	protected internal static bool zeroGroup_1(double[] A, double eps) {
		if ((A[0] < -eps) || (A[0] > eps)) return false;
		if ((A[1] < -eps) || (A[1] > eps)) return false;
		if ((A[2] < -eps) || (A[2] > eps)) return false;
		return true;
	}
	/// <summary>copies coordinates of group 2
	/// </summary>
	protected internal static void copyGroup_2(double[] A, double[] C) {
		copyGroup_1(A, C);
	}
	/// <summary>copies and multiplies (by s) coordinates of group 2
	/// </summary>
	protected internal static void copyMul_2(double[] A, double[] C, double s) {
		copyMul_1(A, C, s);
	}
	/// <summary>copies and divides (by s) coordinates of group 2
	/// </summary>
	protected internal static void copyDiv_2(double[] A, double[] C, double s) {
		copyDiv_1(A, C, s);
	}
	/// <summary>adds coordinates of group 2 from variable A to C
	/// </summary>
	protected internal static void add_2(double[] A, double[] C) {
		add_1(A, C);
	}
	/// <summary>subtracts coordinates of group 2 in variable A from C
	/// </summary>
	protected internal static void sub_2(double[] A, double[] C) {
		sub_1(A, C);
	}
	/// <summary>negate coordinates of group 2 of variable A
	/// </summary>
	protected internal static void neg_2(double[] A, double[] C) {
		neg_1(A, C);
	}
	/// <summary>adds coordinates of group 2 of variables A and B
	/// </summary>
	protected internal static void add2_2_2(double[] A, double[] B, double[] C) {
		add2_1_1(A, B, C);
	}
	/// <summary>subtracts coordinates of group 2 of variables A from B
	/// </summary>
	protected internal static void sub2_2_2(double[] A, double[] B, double[] C) {
		sub2_1_1(A, B, C);
	}
	/// <summary>performs coordinate-wise multiplication of coordinates of group 2 of variables A and B
	/// </summary>
	protected internal static void hp_2_2(double[] A, double[] B, double[] C) {
		hp_1_1(A, B, C);
	}
	/// <summary>performs coordinate-wise division of coordinates of group 2 of variables A and B
	/// (no checks for divide by zero are made)
	/// </summary>
	protected internal static void ihp_2_2(double[] A, double[] B, double[] C) {
		ihp_1_1(A, B, C);
	}
	/// <summary>check for equality up to eps of coordinates of group 2 of variables A and B
	/// </summary>
	protected internal static bool equals_2_2(double[] A, double[] B, double eps) {
		return equals_1_1(A, B, eps);
	}
	/// <summary>checks if coordinates of group 2 of variable A are zero up to eps
	/// </summary>
	protected internal static bool zeroGroup_2(double[] A, double eps) {
		return zeroGroup_1(A, eps);
	}
	/// <summary>copies coordinates of group 3
	/// </summary>
	protected internal static void copyGroup_3(double[] A, double[] C) {
		copyGroup_0(A, C);
	}
	/// <summary>copies and multiplies (by s) coordinates of group 3
	/// </summary>
	protected internal static void copyMul_3(double[] A, double[] C, double s) {
		copyMul_0(A, C, s);
	}
	/// <summary>copies and divides (by s) coordinates of group 3
	/// </summary>
	protected internal static void copyDiv_3(double[] A, double[] C, double s) {
		copyDiv_0(A, C, s);
	}
	/// <summary>adds coordinates of group 3 from variable A to C
	/// </summary>
	protected internal static void add_3(double[] A, double[] C) {
		add_0(A, C);
	}
	/// <summary>subtracts coordinates of group 3 in variable A from C
	/// </summary>
	protected internal static void sub_3(double[] A, double[] C) {
		sub_0(A, C);
	}
	/// <summary>negate coordinates of group 3 of variable A
	/// </summary>
	protected internal static void neg_3(double[] A, double[] C) {
		neg_0(A, C);
	}
	/// <summary>adds coordinates of group 3 of variables A and B
	/// </summary>
	protected internal static void add2_3_3(double[] A, double[] B, double[] C) {
		add2_0_0(A, B, C);
	}
	/// <summary>subtracts coordinates of group 3 of variables A from B
	/// </summary>
	protected internal static void sub2_3_3(double[] A, double[] B, double[] C) {
		sub2_0_0(A, B, C);
	}
	/// <summary>performs coordinate-wise multiplication of coordinates of group 3 of variables A and B
	/// </summary>
	protected internal static void hp_3_3(double[] A, double[] B, double[] C) {
		hp_0_0(A, B, C);
	}
	/// <summary>performs coordinate-wise division of coordinates of group 3 of variables A and B
	/// (no checks for divide by zero are made)
	/// </summary>
	protected internal static void ihp_3_3(double[] A, double[] B, double[] C) {
		ihp_0_0(A, B, C);
	}
	/// <summary>check for equality up to eps of coordinates of group 3 of variables A and B
	/// </summary>
	protected internal static bool equals_3_3(double[] A, double[] B, double eps) {
		return equals_0_0(A, B, eps);
	}
	/// <summary>checks if coordinates of group 3 of variable A are zero up to eps
	/// </summary>
	protected internal static bool zeroGroup_3(double[] A, double eps) {
		return zeroGroup_0(A, eps);
	}
	/// <summary>Computes the partial dual (w.r.t. full space) of a multivector.
	/// </summary>
	protected internal static void dual_default_0_3(double[] A, double[] C) {
		C[0] = -A[0];
	}
	/// <summary>Computes the partial undual (w.r.t. full space) of a multivector.
	/// </summary>
	protected internal static void undual_default_0_3(double[] A, double[] C) {
		C[0] = A[0];
	}
	/// <summary>Computes the partial dual (w.r.t. full space) of a multivector.
	/// </summary>
	protected internal static void dual_default_1_2(double[] A, double[] C) {
		C[0] = -A[2];
		C[1] = A[1];
		C[2] = -A[0];
	}
	/// <summary>Computes the partial undual (w.r.t. full space) of a multivector.
	/// </summary>
	protected internal static void undual_default_1_2(double[] A, double[] C) {
		C[0] = A[2];
		C[1] = -A[1];
		C[2] = A[0];
	}
	/// <summary>Computes the partial dual (w.r.t. full space) of a multivector.
	/// </summary>
	protected internal static void dual_default_2_1(double[] A, double[] C) {
		undual_default_1_2(A, C);
	}
	/// <summary>Computes the partial undual (w.r.t. full space) of a multivector.
	/// </summary>
	protected internal static void undual_default_2_1(double[] A, double[] C) {
		dual_default_1_2(A, C);
	}
	/// <summary>Computes the partial dual (w.r.t. full space) of a multivector.
	/// </summary>
	protected internal static void dual_default_3_0(double[] A, double[] C) {
		undual_default_0_3(A, C);
	}
	/// <summary>Computes the partial undual (w.r.t. full space) of a multivector.
	/// </summary>
	protected internal static void undual_default_3_0(double[] A, double[] C) {
		dual_default_0_3(A, C);
	}
	/// <summary>Computes the partial application of a general outermorphism to a general multivector
	/// </summary>
	protected internal static void applyGomGmv_1_1(om O, double[] A, double[] C) {
		C[0] = (A[0]*O.m_m1[0]+A[1]*O.m_m1[1]+A[2]*O.m_m1[2]);
		C[1] = (A[0]*O.m_m1[3]+A[1]*O.m_m1[4]+A[2]*O.m_m1[5]);
		C[2] = (A[0]*O.m_m1[6]+A[1]*O.m_m1[7]+A[2]*O.m_m1[8]);
	}
	/// <summary>Computes the partial application of a general outermorphism to a general multivector
	/// </summary>
	protected internal static void applyGomGmv_2_2(om O, double[] A, double[] C) {
		C[0] = (A[0]*O.m_m2[0]+A[1]*O.m_m2[1]+A[2]*O.m_m2[2]);
		C[1] = (A[0]*O.m_m2[3]+A[1]*O.m_m2[4]+A[2]*O.m_m2[5]);
		C[2] = (A[0]*O.m_m2[6]+A[1]*O.m_m2[7]+A[2]*O.m_m2[8]);
	}
	/// <summary>Computes the partial application of a general outermorphism to a general multivector
	/// </summary>
	protected internal static void applyGomGmv_3_3(om O, double[] A, double[] C) {
		C[0] = A[0]*O.m_m3[0];
	}
public static bivector _bivector(rotor R)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			R.m_e1_e2, // e1_e2
			R.m_e2_e3, // e2_e3
			R.m_e3_e1 // e3_e1
		);

}
public static vector _vector(oddVersor V)
{
	return new vector(vector.coord_e1_e2_e3,
			V.m_e1, // e1
			V.m_e2, // e2
			V.m_e3 // e3
		);

}
public static trivector _trivector(oddVersor V)
{
	return new trivector(trivector.coord_e1e2e3,
			V.m_e1_e2_e3 // e1_e2_e3
		);

}
/// <summary>Generates a random double in [0.0 1.0) interval using the System.Random class
/// </summary>
public static double genrand() {
	return (double)NextRandomDouble();
}
public static void genrand_seed(int seed) {
	s_randomGenerator = new System.Random(seed);
}

public static void genrand_timeSeed() {
	genrand_seed((int)DateTime.Now.Ticks);
}

/// <summary>Returns mv + mv.
/// </summary>
public static mv add(mv_if a, mv_if b)
{
	double[][] ac = a.to_mv().c();
	double[][] bc = b.to_mv().c();
	double[][] cc = new double[4][];
	
	if (ac[0] != null) {
		cc[0] = new double[1];
		if (bc[0] != null) {
			add2_0_0(ac[0], bc[0], cc[0]);
		}
		else copyGroup_0(ac[0], cc[0]);
	}
	else if (bc[0] != null) {
		cc[0] = new double[1];
		copyGroup_0(bc[0], cc[0]);
	}
	
	if (ac[1] != null) {
		cc[1] = new double[3];
		if (bc[1] != null) {
			add2_1_1(ac[1], bc[1], cc[1]);
		}
		else copyGroup_1(ac[1], cc[1]);
	}
	else if (bc[1] != null) {
		cc[1] = new double[3];
		copyGroup_1(bc[1], cc[1]);
	}
	
	if (ac[2] != null) {
		cc[2] = new double[3];
		if (bc[2] != null) {
			add2_2_2(ac[2], bc[2], cc[2]);
		}
		else copyGroup_2(ac[2], cc[2]);
	}
	else if (bc[2] != null) {
		cc[2] = new double[3];
		copyGroup_2(bc[2], cc[2]);
	}
	
	if (ac[3] != null) {
		cc[3] = new double[1];
		if (bc[3] != null) {
			add2_3_3(ac[3], bc[3], cc[3]);
		}
		else copyGroup_3(ac[3], cc[3]);
	}
	else if (bc[3] != null) {
		cc[3] = new double[1];
		copyGroup_3(bc[3], cc[3]);
	}
	return new mv(cc);
}
/// <summary>Returns vector + vector.
/// </summary>
public static vector add(vector a, vector b)
{
	return new vector(vector.coord_e1_e2_e3,
			(a.m_e1+b.m_e1), // e1
			(a.m_e2+b.m_e2), // e2
			(a.m_e3+b.m_e3) // e3
		);

}
/// <summary>Returns bivector + bivector.
/// </summary>
public static bivector add(bivector a, bivector b)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			(a.m_e1_e2+b.m_e1_e2), // e1_e2
			(a.m_e2_e3+b.m_e2_e3), // e2_e3
			-(-a.m_e3_e1-b.m_e3_e1) // e3_e1
		);

}
/// <summary>Returns vector + trivector.
/// </summary>
public static oddVersor add(vector a, trivector b)
{
	return new oddVersor(oddVersor.coord_e1_e2_e3_e1e2e3,
			a.m_e1, // e1
			a.m_e2, // e2
			a.m_e3, // e3
			b.m_e1_e2_e3 // e1_e2_e3
		);

}
/// <summary>Returns rotor + bivector.
/// </summary>
public static rotor add(rotor a, bivector b)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			a.m_scalar, // scalar
			(a.m_e1_e2+b.m_e1_e2), // e1_e2
			(a.m_e2_e3+b.m_e2_e3), // e2_e3
			-(-a.m_e3_e1-b.m_e3_e1) // e3_e1
		);

}
/// <summary>Returns e1_t + e2_t.
/// </summary>
public static vector add(e1_t a, e2_t b)
{
	return new vector(vector.coord_e1_e2_e3,
			1.0, // e1
			1.0, // e2
			0.0 // e3
		);

}
/// <summary>Returns e1_t + I3_t.
/// </summary>
public static oddVersor add(e1_t a, I3_t b)
{
	return new oddVersor(oddVersor.coord_e1_e2_e3_e1e2e3,
			1.0, // e1
			0.0, // e2
			0.0, // e3
			1.0 // e1_e2_e3
		);

}
/// <summary>Returns mv - mv.
/// </summary>
public static mv subtract(mv_if a, mv_if b)
{
	double[][] ac = a.to_mv().c();
	double[][] bc = b.to_mv().c();
	double[][] cc = new double[4][];
	
	if (ac[0] != null) {
		cc[0] = new double[1];
		if (bc[0] != null) {
			sub2_0_0(ac[0], bc[0], cc[0]);
		}
		else copyGroup_0(ac[0], cc[0]);
	}
	else if (bc[0] != null) {
		cc[0] = new double[1];
		neg_0(bc[0], cc[0]);
	}
	
	if (ac[1] != null) {
		cc[1] = new double[3];
		if (bc[1] != null) {
			sub2_1_1(ac[1], bc[1], cc[1]);
		}
		else copyGroup_1(ac[1], cc[1]);
	}
	else if (bc[1] != null) {
		cc[1] = new double[3];
		neg_1(bc[1], cc[1]);
	}
	
	if (ac[2] != null) {
		cc[2] = new double[3];
		if (bc[2] != null) {
			sub2_2_2(ac[2], bc[2], cc[2]);
		}
		else copyGroup_2(ac[2], cc[2]);
	}
	else if (bc[2] != null) {
		cc[2] = new double[3];
		neg_2(bc[2], cc[2]);
	}
	
	if (ac[3] != null) {
		cc[3] = new double[1];
		if (bc[3] != null) {
			sub2_3_3(ac[3], bc[3], cc[3]);
		}
		else copyGroup_3(ac[3], cc[3]);
	}
	else if (bc[3] != null) {
		cc[3] = new double[1];
		neg_3(bc[3], cc[3]);
	}
	return new mv(cc);
}
/// <summary>Returns vector - vector.
/// </summary>
public static vector subtract(vector a, vector b)
{
	return new vector(vector.coord_e1_e2_e3,
			(a.m_e1-b.m_e1), // e1
			(a.m_e2-b.m_e2), // e2
			(a.m_e3-b.m_e3) // e3
		);

}
/// <summary>Returns bivector - bivector.
/// </summary>
public static bivector subtract(bivector a, bivector b)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			(a.m_e1_e2-b.m_e1_e2), // e1_e2
			(a.m_e2_e3-b.m_e2_e3), // e2_e3
			-(-a.m_e3_e1+b.m_e3_e1) // e3_e1
		);

}
/// <summary>Returns bivector - rotor.
/// </summary>
public static rotor subtract(bivector a, rotor b)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			-b.m_scalar, // scalar
			(a.m_e1_e2-b.m_e1_e2), // e1_e2
			(a.m_e2_e3-b.m_e2_e3), // e2_e3
			-(-a.m_e3_e1+b.m_e3_e1) // e3_e1
		);

}
/// <summary>Returns vector - trivector.
/// </summary>
public static oddVersor subtract(vector a, trivector b)
{
	return new oddVersor(oddVersor.coord_e1_e2_e3_e1e2e3,
			a.m_e1, // e1
			a.m_e2, // e2
			a.m_e3, // e3
			-b.m_e1_e2_e3 // e1_e2_e3
		);

}
/// <summary>Returns a * b * inverse(a) using default metric.
/// </summary>
public static mv applyVersor(mv_if a, mv_if b)
{
	return extractGrade(gp(gp(a, b), versorInverse(a)), b.to_mv().gu());
}
/// <summary>Returns a * b * reverse(a) using default metric. Only gives the correct result when the versor has a positive squared norm.
/// 
/// </summary>
public static mv applyUnitVersor(mv_if a, mv_if b)
{
	return extractGrade(gp(gp(a, b), reverse(a)), b.to_mv().gu());
}
/// <summary>Returns a * b * reverse(a) using default metric. Only gives the correct result when the versor has a positive squared norm.
/// 
/// </summary>
public static mv applyVersorWI(mv_if a, mv_if b, mv_if c)
{
	return extractGrade(gp(gp(a, b), c), b.to_mv().gu());
}
/// <summary>Returns a * b * inverse(a) using default metric.
/// </summary>
public static vector applyVersor(rotor a, vector b)
{
	double _n2_ = (a.m_e1_e2*a.m_e1_e2+a.m_e2_e3*a.m_e2_e3+a.m_e3_e1*a.m_e3_e1+a.m_scalar*a.m_scalar);

	return new vector(vector.coord_e1_e2_e3,
			((-a.m_e1_e2*a.m_e1_e2*b.m_e1+2.0*a.m_e1_e2*a.m_e2_e3*b.m_e3+2.0*a.m_e1_e2*a.m_scalar*b.m_e2+a.m_e2_e3*a.m_e2_e3*b.m_e1+2.0*a.m_e2_e3*a.m_e3_e1*b.m_e2-a.m_e3_e1*a.m_e3_e1*b.m_e1+-2.0*a.m_e3_e1*a.m_scalar*b.m_e3+a.m_scalar*a.m_scalar*b.m_e1))/(_n2_), // e1
			((-a.m_e1_e2*a.m_e1_e2*b.m_e2+2.0*a.m_e1_e2*a.m_e3_e1*b.m_e3+-2.0*a.m_e1_e2*a.m_scalar*b.m_e1-a.m_e2_e3*a.m_e2_e3*b.m_e2+2.0*a.m_e2_e3*a.m_e3_e1*b.m_e1+2.0*a.m_e2_e3*a.m_scalar*b.m_e3+a.m_e3_e1*a.m_e3_e1*b.m_e2+a.m_scalar*a.m_scalar*b.m_e2))/(_n2_), // e2
			((a.m_e1_e2*a.m_e1_e2*b.m_e3+2.0*a.m_e1_e2*a.m_e2_e3*b.m_e1+2.0*a.m_e1_e2*a.m_e3_e1*b.m_e2-a.m_e2_e3*a.m_e2_e3*b.m_e3+-2.0*a.m_e2_e3*a.m_scalar*b.m_e2-a.m_e3_e1*a.m_e3_e1*b.m_e3+2.0*a.m_e3_e1*a.m_scalar*b.m_e1+a.m_scalar*a.m_scalar*b.m_e3))/(_n2_) // e3
		);
}
/// <summary>Returns a * b * reverse(a) using default metric. Only gives the correct result when the versor has a positive squared norm.
/// 
/// </summary>
public static vector applyUnitVersor(rotor a, vector b)
{
	return new vector(vector.coord_e1_e2_e3,
			(-a.m_e1_e2*a.m_e1_e2*b.m_e1+2.0*a.m_e1_e2*a.m_e2_e3*b.m_e3+2.0*a.m_e1_e2*a.m_scalar*b.m_e2+a.m_e2_e3*a.m_e2_e3*b.m_e1+2.0*a.m_e2_e3*a.m_e3_e1*b.m_e2-a.m_e3_e1*a.m_e3_e1*b.m_e1+-2.0*a.m_e3_e1*a.m_scalar*b.m_e3+a.m_scalar*a.m_scalar*b.m_e1), // e1
			(-a.m_e1_e2*a.m_e1_e2*b.m_e2+2.0*a.m_e1_e2*a.m_e3_e1*b.m_e3+-2.0*a.m_e1_e2*a.m_scalar*b.m_e1-a.m_e2_e3*a.m_e2_e3*b.m_e2+2.0*a.m_e2_e3*a.m_e3_e1*b.m_e1+2.0*a.m_e2_e3*a.m_scalar*b.m_e3+a.m_e3_e1*a.m_e3_e1*b.m_e2+a.m_scalar*a.m_scalar*b.m_e2), // e2
			(a.m_e1_e2*a.m_e1_e2*b.m_e3+2.0*a.m_e1_e2*a.m_e2_e3*b.m_e1+2.0*a.m_e1_e2*a.m_e3_e1*b.m_e2-a.m_e2_e3*a.m_e2_e3*b.m_e3+-2.0*a.m_e2_e3*a.m_scalar*b.m_e2-a.m_e3_e1*a.m_e3_e1*b.m_e3+2.0*a.m_e3_e1*a.m_scalar*b.m_e1+a.m_scalar*a.m_scalar*b.m_e3) // e3
		);
}
/// <summary>Returns a * b * reverse(a) using default metric. Only gives the correct result when the versor has a positive squared norm.
/// 
/// </summary>
public static vector applyVersorWI(rotor a, vector b, rotor c)
{
	return new vector(vector.coord_e1_e2_e3,
			(a.m_e1_e2*b.m_e1*c.m_e1_e2+a.m_e1_e2*b.m_e2*c.m_scalar-a.m_e1_e2*b.m_e3*c.m_e2_e3-a.m_e2_e3*b.m_e1*c.m_e2_e3-a.m_e2_e3*b.m_e2*c.m_e3_e1-a.m_e2_e3*b.m_e3*c.m_e1_e2+a.m_e3_e1*b.m_e1*c.m_e3_e1-a.m_e3_e1*b.m_e2*c.m_e2_e3-a.m_e3_e1*b.m_e3*c.m_scalar+a.m_scalar*b.m_e1*c.m_scalar-a.m_scalar*b.m_e2*c.m_e1_e2+a.m_scalar*b.m_e3*c.m_e3_e1), // e1
			(-a.m_e1_e2*b.m_e1*c.m_scalar+a.m_e1_e2*b.m_e2*c.m_e1_e2-a.m_e1_e2*b.m_e3*c.m_e3_e1-a.m_e2_e3*b.m_e1*c.m_e3_e1+a.m_e2_e3*b.m_e2*c.m_e2_e3+a.m_e2_e3*b.m_e3*c.m_scalar-a.m_e3_e1*b.m_e1*c.m_e2_e3-a.m_e3_e1*b.m_e2*c.m_e3_e1-a.m_e3_e1*b.m_e3*c.m_e1_e2+a.m_scalar*b.m_e1*c.m_e1_e2+a.m_scalar*b.m_e2*c.m_scalar-a.m_scalar*b.m_e3*c.m_e2_e3), // e2
			(-a.m_e1_e2*b.m_e1*c.m_e2_e3-a.m_e1_e2*b.m_e2*c.m_e3_e1-a.m_e1_e2*b.m_e3*c.m_e1_e2-a.m_e2_e3*b.m_e1*c.m_e1_e2-a.m_e2_e3*b.m_e2*c.m_scalar+a.m_e2_e3*b.m_e3*c.m_e2_e3+a.m_e3_e1*b.m_e1*c.m_scalar-a.m_e3_e1*b.m_e2*c.m_e1_e2+a.m_e3_e1*b.m_e3*c.m_e3_e1-a.m_scalar*b.m_e1*c.m_e3_e1+a.m_scalar*b.m_e2*c.m_e2_e3+a.m_scalar*b.m_e3*c.m_scalar) // e3
		);
}
/// <summary>Returns a * b * inverse(a) using default metric.
/// </summary>
public static bivector applyVersor(rotor a, bivector b)
{
	double _n2_ = (a.m_e1_e2*a.m_e1_e2+a.m_e2_e3*a.m_e2_e3+a.m_e3_e1*a.m_e3_e1+a.m_scalar*a.m_scalar);

	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			((a.m_e1_e2*a.m_e1_e2*b.m_e1_e2+2.0*a.m_e1_e2*a.m_e2_e3*b.m_e2_e3+2.0*a.m_e1_e2*a.m_e3_e1*b.m_e3_e1-a.m_e2_e3*a.m_e2_e3*b.m_e1_e2+-2.0*a.m_e2_e3*a.m_scalar*b.m_e3_e1-a.m_e3_e1*a.m_e3_e1*b.m_e1_e2+2.0*a.m_e3_e1*a.m_scalar*b.m_e2_e3+a.m_scalar*a.m_scalar*b.m_e1_e2))/(_n2_), // e1_e2
			((-a.m_e1_e2*a.m_e1_e2*b.m_e2_e3+2.0*a.m_e1_e2*a.m_e2_e3*b.m_e1_e2+2.0*a.m_e1_e2*a.m_scalar*b.m_e3_e1+a.m_e2_e3*a.m_e2_e3*b.m_e2_e3+2.0*a.m_e2_e3*a.m_e3_e1*b.m_e3_e1-a.m_e3_e1*a.m_e3_e1*b.m_e2_e3+-2.0*a.m_e3_e1*a.m_scalar*b.m_e1_e2+a.m_scalar*a.m_scalar*b.m_e2_e3))/(_n2_), // e2_e3
			(-(a.m_e1_e2*a.m_e1_e2*b.m_e3_e1+-2.0*a.m_e1_e2*a.m_e3_e1*b.m_e1_e2+2.0*a.m_e1_e2*a.m_scalar*b.m_e2_e3+a.m_e2_e3*a.m_e2_e3*b.m_e3_e1+-2.0*a.m_e2_e3*a.m_e3_e1*b.m_e2_e3+-2.0*a.m_e2_e3*a.m_scalar*b.m_e1_e2-a.m_e3_e1*a.m_e3_e1*b.m_e3_e1-a.m_scalar*a.m_scalar*b.m_e3_e1))/(_n2_) // e3_e1
		);
}
/// <summary>Returns a * b * reverse(a) using default metric. Only gives the correct result when the versor has a positive squared norm.
/// 
/// </summary>
public static bivector applyUnitVersor(rotor a, bivector b)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			(a.m_e1_e2*a.m_e1_e2*b.m_e1_e2+2.0*a.m_e1_e2*a.m_e2_e3*b.m_e2_e3+2.0*a.m_e1_e2*a.m_e3_e1*b.m_e3_e1-a.m_e2_e3*a.m_e2_e3*b.m_e1_e2+-2.0*a.m_e2_e3*a.m_scalar*b.m_e3_e1-a.m_e3_e1*a.m_e3_e1*b.m_e1_e2+2.0*a.m_e3_e1*a.m_scalar*b.m_e2_e3+a.m_scalar*a.m_scalar*b.m_e1_e2), // e1_e2
			(-a.m_e1_e2*a.m_e1_e2*b.m_e2_e3+2.0*a.m_e1_e2*a.m_e2_e3*b.m_e1_e2+2.0*a.m_e1_e2*a.m_scalar*b.m_e3_e1+a.m_e2_e3*a.m_e2_e3*b.m_e2_e3+2.0*a.m_e2_e3*a.m_e3_e1*b.m_e3_e1-a.m_e3_e1*a.m_e3_e1*b.m_e2_e3+-2.0*a.m_e3_e1*a.m_scalar*b.m_e1_e2+a.m_scalar*a.m_scalar*b.m_e2_e3), // e2_e3
			-(a.m_e1_e2*a.m_e1_e2*b.m_e3_e1+-2.0*a.m_e1_e2*a.m_e3_e1*b.m_e1_e2+2.0*a.m_e1_e2*a.m_scalar*b.m_e2_e3+a.m_e2_e3*a.m_e2_e3*b.m_e3_e1+-2.0*a.m_e2_e3*a.m_e3_e1*b.m_e2_e3+-2.0*a.m_e2_e3*a.m_scalar*b.m_e1_e2-a.m_e3_e1*a.m_e3_e1*b.m_e3_e1-a.m_scalar*a.m_scalar*b.m_e3_e1) // e3_e1
		);
}
/// <summary>Returns a * b * reverse(a) using default metric. Only gives the correct result when the versor has a positive squared norm.
/// 
/// </summary>
public static bivector applyVersorWI(rotor a, bivector b, rotor c)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			(-a.m_e1_e2*b.m_e1_e2*c.m_e1_e2-a.m_e1_e2*b.m_e2_e3*c.m_e2_e3-a.m_e1_e2*b.m_e3_e1*c.m_e3_e1+a.m_e2_e3*b.m_e1_e2*c.m_e2_e3-a.m_e2_e3*b.m_e2_e3*c.m_e1_e2-a.m_e2_e3*b.m_e3_e1*c.m_scalar+a.m_e3_e1*b.m_e1_e2*c.m_e3_e1+a.m_e3_e1*b.m_e2_e3*c.m_scalar-a.m_e3_e1*b.m_e3_e1*c.m_e1_e2+a.m_scalar*b.m_e1_e2*c.m_scalar-a.m_scalar*b.m_e2_e3*c.m_e3_e1+a.m_scalar*b.m_e3_e1*c.m_e2_e3), // e1_e2
			(-a.m_e1_e2*b.m_e1_e2*c.m_e2_e3+a.m_e1_e2*b.m_e2_e3*c.m_e1_e2+a.m_e1_e2*b.m_e3_e1*c.m_scalar-a.m_e2_e3*b.m_e1_e2*c.m_e1_e2-a.m_e2_e3*b.m_e2_e3*c.m_e2_e3-a.m_e2_e3*b.m_e3_e1*c.m_e3_e1-a.m_e3_e1*b.m_e1_e2*c.m_scalar+a.m_e3_e1*b.m_e2_e3*c.m_e3_e1-a.m_e3_e1*b.m_e3_e1*c.m_e2_e3+a.m_scalar*b.m_e1_e2*c.m_e3_e1+a.m_scalar*b.m_e2_e3*c.m_scalar-a.m_scalar*b.m_e3_e1*c.m_e1_e2), // e2_e3
			-(a.m_e1_e2*b.m_e1_e2*c.m_e3_e1+a.m_e1_e2*b.m_e2_e3*c.m_scalar-a.m_e1_e2*b.m_e3_e1*c.m_e1_e2-a.m_e2_e3*b.m_e1_e2*c.m_scalar+a.m_e2_e3*b.m_e2_e3*c.m_e3_e1-a.m_e2_e3*b.m_e3_e1*c.m_e2_e3+a.m_e3_e1*b.m_e1_e2*c.m_e1_e2+a.m_e3_e1*b.m_e2_e3*c.m_e2_e3+a.m_e3_e1*b.m_e3_e1*c.m_e3_e1+a.m_scalar*b.m_e1_e2*c.m_e2_e3-a.m_scalar*b.m_e2_e3*c.m_e1_e2-a.m_scalar*b.m_e3_e1*c.m_scalar) // e3_e1
		);
}
/// <summary>Returns a * b * inverse(a) using default metric.
/// </summary>
public static trivector applyVersor(rotor a, trivector b)
{
	double _n2_ = (a.m_e1_e2*a.m_e1_e2+a.m_e2_e3*a.m_e2_e3+a.m_e3_e1*a.m_e3_e1+a.m_scalar*a.m_scalar);

	return new trivector(trivector.coord_e1e2e3,
			((a.m_e1_e2*a.m_e1_e2*b.m_e1_e2_e3+a.m_e2_e3*a.m_e2_e3*b.m_e1_e2_e3+a.m_e3_e1*a.m_e3_e1*b.m_e1_e2_e3+a.m_scalar*a.m_scalar*b.m_e1_e2_e3))/(_n2_) // e1_e2_e3
		);
}
/// <summary>Returns a * b * reverse(a) using default metric. Only gives the correct result when the versor has a positive squared norm.
/// 
/// </summary>
public static trivector applyUnitVersor(rotor a, trivector b)
{
	return new trivector(trivector.coord_e1e2e3,
			(a.m_e1_e2*a.m_e1_e2*b.m_e1_e2_e3+a.m_e2_e3*a.m_e2_e3*b.m_e1_e2_e3+a.m_e3_e1*a.m_e3_e1*b.m_e1_e2_e3+a.m_scalar*a.m_scalar*b.m_e1_e2_e3) // e1_e2_e3
		);
}
/// <summary>Returns a * b * reverse(a) using default metric. Only gives the correct result when the versor has a positive squared norm.
/// 
/// </summary>
public static trivector applyVersorWI(rotor a, trivector b, rotor c)
{
	return new trivector(trivector.coord_e1e2e3,
			(-a.m_e1_e2*b.m_e1_e2_e3*c.m_e1_e2-a.m_e2_e3*b.m_e1_e2_e3*c.m_e2_e3-a.m_e3_e1*b.m_e1_e2_e3*c.m_e3_e1+a.m_scalar*b.m_e1_e2_e3*c.m_scalar) // e1_e2_e3
		);
}
/// <summary>Returns a * b * inverse(a) using default metric.
/// </summary>
public static vector applyVersor(rotor a, e1_t b)
{
	double _n2_ = (a.m_e1_e2*a.m_e1_e2+a.m_e2_e3*a.m_e2_e3+a.m_e3_e1*a.m_e3_e1+a.m_scalar*a.m_scalar);

	return new vector(vector.coord_e1_e2_e3,
			((-a.m_e1_e2*a.m_e1_e2+a.m_e2_e3*a.m_e2_e3-a.m_e3_e1*a.m_e3_e1+a.m_scalar*a.m_scalar))/(_n2_), // e1
			((-2.0*a.m_e1_e2*a.m_scalar+2.0*a.m_e2_e3*a.m_e3_e1))/(_n2_), // e2
			((2.0*a.m_e1_e2*a.m_e2_e3+2.0*a.m_e3_e1*a.m_scalar))/(_n2_) // e3
		);
}
/// <summary>Returns a * b * reverse(a) using default metric. Only gives the correct result when the versor has a positive squared norm.
/// 
/// </summary>
public static vector applyUnitVersor(rotor a, e2_t b)
{
	return new vector(vector.coord_e1_e2_e3,
			(2.0*a.m_e1_e2*a.m_scalar+2.0*a.m_e2_e3*a.m_e3_e1), // e1
			(-a.m_e1_e2*a.m_e1_e2-a.m_e2_e3*a.m_e2_e3+a.m_e3_e1*a.m_e3_e1+a.m_scalar*a.m_scalar), // e2
			(2.0*a.m_e1_e2*a.m_e3_e1+-2.0*a.m_e2_e3*a.m_scalar) // e3
		);
}
/// <summary>Returns a * b * inverse(a) using default metric.
/// </summary>
public static trivector applyVersor(rotor a, I3_t b)
{
	double _n2_ = (a.m_e1_e2*a.m_e1_e2+a.m_e2_e3*a.m_e2_e3+a.m_e3_e1*a.m_e3_e1+a.m_scalar*a.m_scalar);

	return new trivector(trivector.coord_e1e2e3,
			((a.m_e1_e2*a.m_e1_e2+a.m_e2_e3*a.m_e2_e3+a.m_e3_e1*a.m_e3_e1+a.m_scalar*a.m_scalar))/(_n2_) // e1_e2_e3
		);
}
/// <summary>Returns a * b * reverse(a) using default metric. Only gives the correct result when the versor has a positive squared norm.
/// 
/// </summary>
public static trivector applyUnitVersor(rotor a, I3_t b)
{
	return new trivector(trivector.coord_e1e2e3,
			(a.m_e1_e2*a.m_e1_e2+a.m_e2_e3*a.m_e2_e3+a.m_e3_e1*a.m_e3_e1+a.m_scalar*a.m_scalar) // e1_e2_e3
		);
}
/// <summary>Returns a * b * reverse(a) using default metric. Only gives the correct result when the versor has a positive squared norm.
/// 
/// </summary>
public static trivector applyVersorWI(rotor a, I3_t b, rotor c)
{
	return new trivector(trivector.coord_e1e2e3,
			(-a.m_e1_e2*c.m_e1_e2-a.m_e2_e3*c.m_e2_e3-a.m_e3_e1*c.m_e3_e1+a.m_scalar*c.m_scalar) // e1_e2_e3
		);
}
/// <summary>Returns om * mv.
/// </summary>
public static mv applyOM(om a, mv_if b)
{
	double[][] bc = b.to_mv().c();
	double[][] cc = new double[4][];
	if (bc[0] != null) {
	}
	
	if (bc[1] != null) {
		if (cc[1] == null) cc[1] = new double[3];
		applyGomGmv_1_1(a, bc[1], cc[1]);
	}
	
	if (bc[2] != null) {
		if (cc[2] == null) cc[2] = new double[3];
		applyGomGmv_2_2(a, bc[2], cc[2]);
	}
	
	if (bc[3] != null) {
		if (cc[3] == null) cc[3] = new double[1];
		applyGomGmv_3_3(a, bc[3], cc[3]);
	}
	
	return new mv(cc);
}
/// <summary>Returns om * vector.
/// </summary>
public static vector applyOM(om a, vector b)
{
	return new vector(vector.coord_e1_e2_e3,
			(a.m_m1[0]*b.m_e1+a.m_m1[1]*b.m_e2+a.m_m1[2]*b.m_e3), // e1
			(a.m_m1[3]*b.m_e1+a.m_m1[4]*b.m_e2+a.m_m1[5]*b.m_e3), // e2
			(a.m_m1[6]*b.m_e1+a.m_m1[7]*b.m_e2+a.m_m1[8]*b.m_e3) // e3
		);

}
/// <summary>Returns om * bivector.
/// </summary>
public static bivector applyOM(om a, bivector b)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			(a.m_m2[0]*b.m_e1_e2-a.m_m2[1]*b.m_e3_e1+a.m_m2[2]*b.m_e2_e3), // e1_e2
			(a.m_m2[6]*b.m_e1_e2-a.m_m2[7]*b.m_e3_e1+a.m_m2[8]*b.m_e2_e3), // e2_e3
			-(a.m_m2[3]*b.m_e1_e2-a.m_m2[4]*b.m_e3_e1+a.m_m2[5]*b.m_e2_e3) // e3_e1
		);

}
/// <summary>Returns om * trivector.
/// </summary>
public static trivector applyOM(om a, trivector b)
{
	return new trivector(trivector.coord_e1e2e3,
			a.m_m3[0]*b.m_e1_e2_e3 // e1_e2_e3
		);

}
/// <summary>Returns grade1OM * vector.
/// </summary>
public static vector applyOM(grade1OM a, vector b)
{
	return new vector(vector.coord_e1_e2_e3,
			(a.m_m1[0]*b.m_e1+a.m_m1[1]*b.m_e2+a.m_m1[2]*b.m_e3), // e1
			(a.m_m1[3]*b.m_e1+a.m_m1[4]*b.m_e2+a.m_m1[5]*b.m_e3), // e2
			(a.m_m1[6]*b.m_e1+a.m_m1[7]*b.m_e2+a.m_m1[8]*b.m_e3) // e3
		);

}
/// <summary>Returns grade2OM * bivector.
/// </summary>
public static bivector applyOM(grade2OM a, bivector b)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			(a.m_m2[0]*b.m_e1_e2+a.m_m2[1]*b.m_e2_e3+a.m_m2[2]*b.m_e3_e1), // e1_e2
			(a.m_m2[3]*b.m_e1_e2+a.m_m2[4]*b.m_e2_e3+a.m_m2[5]*b.m_e3_e1), // e2_e3
			-(-a.m_m2[6]*b.m_e1_e2-a.m_m2[7]*b.m_e2_e3-a.m_m2[8]*b.m_e3_e1) // e3_e1
		);

}
/// <summary>Returns grade3OM * trivector.
/// </summary>
public static trivector applyOM(grade3OM a, trivector b)
{
	return new trivector(trivector.coord_e1e2e3,
			a.m_m3[0]*b.m_e1_e2_e3 // e1_e2_e3
		);

}
/// <summary>Returns grade1OM * e1_t.
/// </summary>
public static vector applyOM(grade1OM a, e1_t b)
{
	return new vector(vector.coord_e1_e2_e3,
			a.m_m1[0], // e1
			a.m_m1[3], // e2
			a.m_m1[6] // e3
		);

}
/// <summary>Returns grade1OM * e2_t.
/// </summary>
public static vector applyOM(grade1OM a, e2_t b)
{
	return new vector(vector.coord_e1_e2_e3,
			a.m_m1[1], // e1
			a.m_m1[4], // e2
			a.m_m1[7] // e3
		);

}
/// <summary>Returns grade1OM * e3_t.
/// </summary>
public static vector applyOM(grade1OM a, e3_t b)
{
	return new vector(vector.coord_e1_e2_e3,
			a.m_m1[2], // e1
			a.m_m1[5], // e2
			a.m_m1[8] // e3
		);

}
/// <summary>Returns grade3OM * I3_t.
/// </summary>
public static trivector applyOM(grade3OM a, I3_t b)
{
	return new trivector(trivector.coord_e1e2e3,
			a.m_m3[0] // e1_e2_e3
		);

}
/// <summary>Returns a / b
/// </summary>
public static mv div(mv_if a, double b)
{
	double[][] ac = a.to_mv().c();
	double[][] cc = new double[4][];
	
	if (ac[0] != null) {
		cc[0] = new double[1];
		copyDiv_0(ac[0], cc[0], b);
	}
	
	if (ac[1] != null) {
		cc[1] = new double[3];
		copyDiv_1(ac[1], cc[1], b);
	}
	
	if (ac[2] != null) {
		cc[2] = new double[3];
		copyDiv_2(ac[2], cc[2], b);
	}
	
	if (ac[3] != null) {
		cc[3] = new double[1];
		copyDiv_3(ac[3], cc[3], b);
	}
	return new mv(cc);
}
/// <summary>Returns a / b
/// </summary>
public static vector div(vector a, double b)
{
	return new vector(vector.coord_e1_e2_e3,
			a.m_e1/((b)), // e1
			a.m_e2/((b)), // e2
			a.m_e3/((b)) // e3
		);
}
/// <summary>Returns a / b
/// </summary>
public static bivector div(bivector a, double b)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			a.m_e1_e2/((b)), // e1_e2
			a.m_e2_e3/((b)), // e2_e3
			a.m_e3_e1/((b)) // e3_e1
		);
}
/// <summary>Returns a / b
/// </summary>
public static trivector div(trivector a, double b)
{
	return new trivector(trivector.coord_e1e2e3,
			a.m_e1_e2_e3/((b)) // e1_e2_e3
		);
}
/// <summary>Returns a / b
/// </summary>
public static rotor div(rotor a, double b)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			a.m_scalar/((b)), // scalar
			a.m_e1_e2/((b)), // e1_e2
			a.m_e2_e3/((b)), // e2_e3
			a.m_e3_e1/((b)) // e3_e1
		);
}
/// <summary>Returns a / b
/// </summary>
public static vector div(e1_t a, double b)
{
	return new vector(vector.coord_e1_e2_e3,
			1.0 / (b), // e1
			0.0, // e2
			0.0 // e3
		);
}
/// <summary>Returns a / b
/// </summary>
public static trivector div(I3_t a, double b)
{
	return new trivector(trivector.coord_e1e2e3,
			1.0 / (b) // e1_e2_e3
		);
}
/// <summary>Returns dual of mv using default metric.
/// </summary>
public static mv dual(mv_if a)
{
	double[][] ac = a.to_mv().c();
	double[][] cc = new double[4][];
	if (ac[0] != null) {
		if (cc[3] == null) cc[3] = new double[1];
		dual_default_0_3(ac[0], cc[3]);
	}
	
	if (ac[1] != null) {
		if (cc[2] == null) cc[2] = new double[3];
		dual_default_1_2(ac[1], cc[2]);
	}
	
	if (ac[2] != null) {
		if (cc[1] == null) cc[1] = new double[3];
		dual_default_2_1(ac[2], cc[1]);
	}
	
	if (ac[3] != null) {
		if (cc[0] == null) cc[0] = new double[1];
		dual_default_3_0(ac[3], cc[0]);
	}
	
	return new mv(cc);
}
/// <summary>Returns undual of mv using default metric.
/// </summary>
public static mv undual(mv_if a)
{
	double[][] ac = a.to_mv().c();
	double[][] cc = new double[4][];
	if (ac[0] != null) {
		if (cc[3] == null) cc[3] = new double[1];
		undual_default_0_3(ac[0], cc[3]);
	}
	
	if (ac[1] != null) {
		if (cc[2] == null) cc[2] = new double[3];
		undual_default_1_2(ac[1], cc[2]);
	}
	
	if (ac[2] != null) {
		if (cc[1] == null) cc[1] = new double[3];
		undual_default_2_1(ac[2], cc[1]);
	}
	
	if (ac[3] != null) {
		if (cc[0] == null) cc[0] = new double[1];
		undual_default_3_0(ac[3], cc[0]);
	}
	
	return new mv(cc);
}
/// <summary>Returns dual of double using default metric.
/// </summary>
public static trivector dual(double a)
{
	return new trivector(trivector.coord_e1e2e3,
			-a // e1_e2_e3
		);

}
/// <summary>Returns undual of double using default metric.
/// </summary>
public static trivector undual(double a)
{
	return new trivector(trivector.coord_e1e2e3,
			a // e1_e2_e3
		);

}
/// <summary>Returns dual of vector using default metric.
/// </summary>
public static bivector dual(vector a)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			-a.m_e3, // e1_e2
			-a.m_e1, // e2_e3
			-a.m_e2 // e3_e1
		);

}
/// <summary>Returns undual of vector using default metric.
/// </summary>
public static bivector undual(vector a)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			a.m_e3, // e1_e2
			a.m_e1, // e2_e3
			a.m_e2 // e3_e1
		);

}
/// <summary>Returns dual of bivector using default metric.
/// </summary>
public static vector dual(bivector a)
{
	return new vector(vector.coord_e1_e2_e3,
			a.m_e2_e3, // e1
			a.m_e3_e1, // e2
			a.m_e1_e2 // e3
		);

}
/// <summary>Returns undual of bivector using default metric.
/// </summary>
public static vector undual(bivector a)
{
	return new vector(vector.coord_e1_e2_e3,
			-a.m_e2_e3, // e1
			-a.m_e3_e1, // e2
			-a.m_e1_e2 // e3
		);

}
/// <summary>Returns dual of rotor using default metric.
/// </summary>
public static oddVersor dual(rotor a)
{
	return new oddVersor(oddVersor.coord_e1_e2_e3_e1e2e3,
			a.m_e2_e3, // e1
			a.m_e3_e1, // e2
			a.m_e1_e2, // e3
			-a.m_scalar // e1_e2_e3
		);

}
/// <summary>Returns undual of rotor using default metric.
/// </summary>
public static oddVersor undual(rotor a)
{
	return new oddVersor(oddVersor.coord_e1_e2_e3_e1e2e3,
			-a.m_e2_e3, // e1
			-a.m_e3_e1, // e2
			-a.m_e1_e2, // e3
			a.m_scalar // e1_e2_e3
		);

}
/// <summary>Returns dual of oddVersor using default metric.
/// </summary>
public static rotor dual(oddVersor a)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			a.m_e1_e2_e3, // scalar
			-a.m_e3, // e1_e2
			-a.m_e1, // e2_e3
			-a.m_e2 // e3_e1
		);

}
/// <summary>Returns undual of oddVersor using default metric.
/// </summary>
public static rotor undual(oddVersor a)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			-a.m_e1_e2_e3, // scalar
			a.m_e3, // e1_e2
			a.m_e1, // e2_e3
			a.m_e2 // e3_e1
		);

}
/// <summary>Returns dual of trivector using default metric.
/// </summary>
public static double dual(trivector a)
{
	return a.m_e1_e2_e3;

}
/// <summary>Returns undual of trivector using default metric.
/// </summary>
public static double undual(trivector a)
{
	return -a.m_e1_e2_e3;

}
/// <summary>Returns dual of e1_t using default metric.
/// </summary>
public static bivector dual(e1_t a)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			0.0, // e1_e2
			-1.0, // e2_e3
			0.0 // e3_e1
		);

}
/// <summary>Returns undual of e2_t using default metric.
/// </summary>
public static bivector undual(e2_t a)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			0.0, // e1_e2
			0.0, // e2_e3
			1.0 // e3_e1
		);

}
/// <summary>Returns dual of I3_t using default metric.
/// </summary>
public static double dual(I3_t a)
{
	return 1.0;

}
/// <summary>Returns undual of I3_t using default metric.
/// </summary>
public static double undual(I3_t a)
{
	return -1.0;

}
/// <summary>Returns whether input multivectors are equal up to an epsilon c.
/// </summary>
public static bool equals(mv_if a, mv_if b, double c)
{
	double[][] ac = a.to_mv().c();
	double[][] bc = b.to_mv().c();
	
	if (ac[0] != null) {
		if (bc[0] != null) {
			if (!equals_0_0(ac[0], bc[0], c)) return false;
		}
		else if (!zeroGroup_0(ac[0], c)) return false;
	}
		else if (bc[0] != null) {
		if (!zeroGroup_0(bc[0], c)) return false;
	}
	
	if (ac[1] != null) {
		if (bc[1] != null) {
			if (!equals_1_1(ac[1], bc[1], c)) return false;
		}
		else if (!zeroGroup_1(ac[1], c)) return false;
	}
		else if (bc[1] != null) {
		if (!zeroGroup_1(bc[1], c)) return false;
	}
	
	if (ac[2] != null) {
		if (bc[2] != null) {
			if (!equals_2_2(ac[2], bc[2], c)) return false;
		}
		else if (!zeroGroup_2(ac[2], c)) return false;
	}
		else if (bc[2] != null) {
		if (!zeroGroup_2(bc[2], c)) return false;
	}
	
	if (ac[3] != null) {
		if (bc[3] != null) {
			if (!equals_3_3(ac[3], bc[3], c)) return false;
		}
		else if (!zeroGroup_3(ac[3], c)) return false;
	}
		else if (bc[3] != null) {
		if (!zeroGroup_3(bc[3], c)) return false;
	}
	return true;
}
/// <summary>Returns whether input multivectors are equal up to an epsilon c.
/// </summary>
public static bool equals(vector a, vector b, double c)
{
	double d;
	d = a.m_e1 - b.m_e1; if ((d < -c) || (d > c)) return false; /* e1 */
	d = a.m_e2 - b.m_e2; if ((d < -c) || (d > c)) return false; /* e2 */
	d = a.m_e3 - b.m_e3; if ((d < -c) || (d > c)) return false; /* e3 */
	return true;
}
/// <summary>Returns whether input multivectors are equal up to an epsilon c.
/// </summary>
public static bool equals(bivector a, bivector b, double c)
{
	double d;
	d = a.m_e1_e2 - b.m_e1_e2; if ((d < -c) || (d > c)) return false; /* e1^e2 */
	d = -a.m_e3_e1 - -b.m_e3_e1; if ((d < -c) || (d > c)) return false; /* e1^e3 */
	d = a.m_e2_e3 - b.m_e2_e3; if ((d < -c) || (d > c)) return false; /* e2^e3 */
	return true;
}
/// <summary>Returns whether input multivectors are equal up to an epsilon c.
/// </summary>
public static bool equals(rotor a, rotor b, double c)
{
	double d;
	d = a.m_scalar - b.m_scalar; if ((d < -c) || (d > c)) return false; /* 1 */
	d = a.m_e1_e2 - b.m_e1_e2; if ((d < -c) || (d > c)) return false; /* e1^e2 */
	d = -a.m_e3_e1 - -b.m_e3_e1; if ((d < -c) || (d > c)) return false; /* e1^e3 */
	d = a.m_e2_e3 - b.m_e2_e3; if ((d < -c) || (d > c)) return false; /* e2^e3 */
	return true;
}
/// <summary>Returns whether input multivectors are equal up to an epsilon c.
/// </summary>
public static bool equals(bivector a, rotor b, double c)
{
	double d;
	if ((b.m_scalar < -c) || (b.m_scalar > c)) return false; /* 1 */
	d = a.m_e1_e2 - b.m_e1_e2; if ((d < -c) || (d > c)) return false; /* e1^e2 */
	d = -a.m_e3_e1 - -b.m_e3_e1; if ((d < -c) || (d > c)) return false; /* e1^e3 */
	d = a.m_e2_e3 - b.m_e2_e3; if ((d < -c) || (d > c)) return false; /* e2^e3 */
	return true;
}
/// <summary>Returns whether input multivectors are equal up to an epsilon c.
/// </summary>
public static bool equals(trivector a, trivector b, double c)
{
	double d;
	d = a.m_e1_e2_e3 - b.m_e1_e2_e3; if ((d < -c) || (d > c)) return false; /* e1^e2^e3 */
	return true;
}
/// <summary>Returns whether input multivectors are equal up to an epsilon c.
/// </summary>
public static bool equals(rotor a, bivector b, double c)
{
	double d;
	if ((a.m_scalar < -c) || (a.m_scalar > c)) return false; /* 1 */
	d = a.m_e1_e2 - b.m_e1_e2; if ((d < -c) || (d > c)) return false; /* e1^e2 */
	d = -a.m_e3_e1 - -b.m_e3_e1; if ((d < -c) || (d > c)) return false; /* e1^e3 */
	d = a.m_e2_e3 - b.m_e2_e3; if ((d < -c) || (d > c)) return false; /* e2^e3 */
	return true;
}
/// <summary>Returns whether input multivectors are equal up to an epsilon c.
/// </summary>
public static bool equals(e1_t a, e1_t b, double c)
{
	double d;
	d = 1.0 - 1.0; if ((d < -c) || (d > c)) return false; /* e1 */
	return true;
}
/// <summary>Returns whether input multivectors are equal up to an epsilon c.
/// </summary>
public static bool equals(e2_t a, I3_t b, double c)
{
	if ((1.0 < -c) || (1.0 > c)) return false; /* e2 */
	if ((1.0 < -c) || (1.0 > c)) return false; /* e1^e2^e3 */
	return true;
}
/// <summary>Returns grade groupBitmap of  mv.
/// </summary>
public static mv extractGrade(mv_if a, GroupBitmap groupBitmap)
{
	GroupBitmap gu = a.to_mv().gu() & groupBitmap;
	double[][] ac = a.to_mv().c();
	double[][] cc = new double[4][];
	
	if ((gu & GroupBitmap.GROUP_0) != 0) {
		cc[0] = new double[1];
		copyGroup_0(ac[0], cc[0]);
	}
	
	if ((gu & GroupBitmap.GROUP_1) != 0) {
		cc[1] = new double[3];
		copyGroup_1(ac[1], cc[1]);
	}
	
	if ((gu & GroupBitmap.GROUP_2) != 0) {
		cc[2] = new double[3];
		copyGroup_2(ac[2], cc[2]);
	}
	
	if ((gu & GroupBitmap.GROUP_3) != 0) {
		cc[3] = new double[1];
		copyGroup_3(ac[3], cc[3]);
	}
	return new mv(cc);
}
/// <summary>Returns grade 2 of  mv.
/// </summary>
public static mv extractGrade2(mv_if a)
{
	return extractGrade(a, (GroupBitmap)4);
}
/// <summary>Returns grade 0 of  rotor.
/// </summary>
public static double extractGrade0(rotor a)
{
	return a.m_scalar;
}
/// <summary>Returns grade 2 of  rotor.
/// </summary>
public static bivector extractGrade2(rotor a)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			a.m_e1_e2, // e1_e2
			a.m_e2_e3, // e2_e3
			a.m_e3_e1 // e3_e1
		);
}
/// <summary>Returns grade 0 of  oddVersor.
/// </summary>
public static double extractGrade0(oddVersor a)
{
	return 0.0;
}
/// <summary>Returns grade 1 of  oddVersor.
/// </summary>
public static vector extractGrade1(oddVersor a)
{
	return new vector(vector.coord_e1_e2_e3,
			a.m_e1, // e1
			a.m_e2, // e2
			a.m_e3 // e3
		);
}
/// <summary>Returns grade 2 of  oddVersor.
/// </summary>
public static double extractGrade2(oddVersor a)
{
	return 0.0;
}
/// <summary>Returns grade 3 of  oddVersor.
/// </summary>
public static trivector extractGrade3(oddVersor a)
{
	return new trivector(trivector.coord_e1e2e3,
			a.m_e1_e2_e3 // e1_e2_e3
		);
}
/// <summary>Returns grade 0 of  e1_t.
/// </summary>
public static double extractGrade0(e1_t a)
{
	return 0.0;
}
/// <summary>Returns grade 1 of  e2_t.
/// </summary>
public static e2_t extractGrade1(e2_t a)
{
	return new e2_t(		);
}
/// <summary>Returns grade 2 of  e3_t.
/// </summary>
public static double extractGrade2(e3_t a)
{
	return 0.0;
}
/// <summary>Returns grade 3 of  e1_t.
/// </summary>
public static double extractGrade3(e1_t a)
{
	return 0.0;
}
/// <summary>Returns grade 0 of  I3_t.
/// </summary>
public static double extractGrade0(I3_t a)
{
	return 0.0;
}
/// <summary>Returns grade 1 of  I3_t.
/// </summary>
public static double extractGrade1(I3_t a)
{
	return 0.0;
}
/// <summary>Returns grade 2 of  I3_t.
/// </summary>
public static double extractGrade2(I3_t a)
{
	return 0.0;
}
/// <summary>Returns grade 3 of  I3_t.
/// </summary>
public static I3_t extractGrade3(I3_t a)
{
	return new I3_t(		);
}
/// <summary>Returns geometric product of mv and mv.
/// </summary>
public static mv gp(mv_if a, mv_if b)
{
	double[][] ac = a.to_mv().c();
	double[][] bc = b.to_mv().c();
	double[][] cc = new double[4][];
	if (ac[0] != null) {
		if (bc[0] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_0_0_0(ac[0], bc[0], cc[0]);
		}
		if (bc[1] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_0_1_1(ac[0], bc[1], cc[1]);
		}
		if (bc[2] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_0_2_2(ac[0], bc[2], cc[2]);
		}
		if (bc[3] != null) {
			if (cc[3] == null) cc[3] = new double[1];
			gp_default_0_3_3(ac[0], bc[3], cc[3]);
		}
	}
	if (ac[1] != null) {
		if (bc[0] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_1_0_1(ac[1], bc[0], cc[1]);
		}
		if (bc[1] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_1_1_0(ac[1], bc[1], cc[0]);
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_1_1_2(ac[1], bc[1], cc[2]);
		}
		if (bc[2] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_1_2_1(ac[1], bc[2], cc[1]);
			if (cc[3] == null) cc[3] = new double[1];
			gp_default_1_2_3(ac[1], bc[2], cc[3]);
		}
		if (bc[3] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_1_3_2(ac[1], bc[3], cc[2]);
		}
	}
	if (ac[2] != null) {
		if (bc[0] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_2_0_2(ac[2], bc[0], cc[2]);
		}
		if (bc[1] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_2_1_1(ac[2], bc[1], cc[1]);
			if (cc[3] == null) cc[3] = new double[1];
			gp_default_2_1_3(ac[2], bc[1], cc[3]);
		}
		if (bc[2] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_2_2_0(ac[2], bc[2], cc[0]);
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_2_2_2(ac[2], bc[2], cc[2]);
		}
		if (bc[3] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_2_3_1(ac[2], bc[3], cc[1]);
		}
	}
	if (ac[3] != null) {
		if (bc[0] != null) {
			if (cc[3] == null) cc[3] = new double[1];
			gp_default_3_0_3(ac[3], bc[0], cc[3]);
		}
		if (bc[1] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_3_1_2(ac[3], bc[1], cc[2]);
		}
		if (bc[2] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_3_2_1(ac[3], bc[2], cc[1]);
		}
		if (bc[3] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_3_3_0(ac[3], bc[3], cc[0]);
		}
	}
	return new mv(cc);
}
/// <summary>Returns geometric product of vector and vector.
/// </summary>
public static rotor gp(vector a, vector b)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			(a.m_e1*b.m_e1+a.m_e2*b.m_e2+a.m_e3*b.m_e3), // scalar
			(a.m_e1*b.m_e2-a.m_e2*b.m_e1), // e1_e2
			(a.m_e2*b.m_e3-a.m_e3*b.m_e2), // e2_e3
			-(a.m_e1*b.m_e3-a.m_e3*b.m_e1) // e3_e1
		);

}
/// <summary>Returns geometric product of rotor and vector.
/// </summary>
public static oddVersor gp(rotor a, vector b)
{
	return new oddVersor(oddVersor.coord_e1_e2_e3_e1e2e3,
			(a.m_e1_e2*b.m_e2-a.m_e3_e1*b.m_e3+a.m_scalar*b.m_e1), // e1
			(-a.m_e1_e2*b.m_e1+a.m_e2_e3*b.m_e3+a.m_scalar*b.m_e2), // e2
			(-a.m_e2_e3*b.m_e2+a.m_e3_e1*b.m_e1+a.m_scalar*b.m_e3), // e3
			(a.m_e1_e2*b.m_e3+a.m_e2_e3*b.m_e1+a.m_e3_e1*b.m_e2) // e1_e2_e3
		);

}
/// <summary>Returns geometric product of vector and rotor.
/// </summary>
public static oddVersor gp(vector a, rotor b)
{
	return new oddVersor(oddVersor.coord_e1_e2_e3_e1e2e3,
			(a.m_e1*b.m_scalar-a.m_e2*b.m_e1_e2+a.m_e3*b.m_e3_e1), // e1
			(a.m_e1*b.m_e1_e2+a.m_e2*b.m_scalar-a.m_e3*b.m_e2_e3), // e2
			(-a.m_e1*b.m_e3_e1+a.m_e2*b.m_e2_e3+a.m_e3*b.m_scalar), // e3
			(a.m_e1*b.m_e2_e3+a.m_e2*b.m_e3_e1+a.m_e3*b.m_e1_e2) // e1_e2_e3
		);

}
/// <summary>Returns geometric product of rotor and rotor.
/// </summary>
public static rotor gp(rotor a, rotor b)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			(-a.m_e1_e2*b.m_e1_e2-a.m_e2_e3*b.m_e2_e3-a.m_e3_e1*b.m_e3_e1+a.m_scalar*b.m_scalar), // scalar
			(a.m_e1_e2*b.m_scalar-a.m_e2_e3*b.m_e3_e1+a.m_e3_e1*b.m_e2_e3+a.m_scalar*b.m_e1_e2), // e1_e2
			(a.m_e1_e2*b.m_e3_e1+a.m_e2_e3*b.m_scalar-a.m_e3_e1*b.m_e1_e2+a.m_scalar*b.m_e2_e3), // e2_e3
			-(a.m_e1_e2*b.m_e2_e3-a.m_e2_e3*b.m_e1_e2-a.m_e3_e1*b.m_scalar-a.m_scalar*b.m_e3_e1) // e3_e1
		);

}
/// <summary>Returns geometric product of bivector and bivector.
/// </summary>
public static rotor gp(bivector a, bivector b)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			(-a.m_e1_e2*b.m_e1_e2-a.m_e2_e3*b.m_e2_e3-a.m_e3_e1*b.m_e3_e1), // scalar
			(-a.m_e2_e3*b.m_e3_e1+a.m_e3_e1*b.m_e2_e3), // e1_e2
			(a.m_e1_e2*b.m_e3_e1-a.m_e3_e1*b.m_e1_e2), // e2_e3
			-(a.m_e1_e2*b.m_e2_e3-a.m_e2_e3*b.m_e1_e2) // e3_e1
		);

}
/// <summary>Returns geometric product of e1_t and rotor.
/// </summary>
public static oddVersor gp(e1_t a, rotor b)
{
	return new oddVersor(oddVersor.coord_e1_e2_e3_e1e2e3,
			b.m_scalar, // e1
			b.m_e1_e2, // e2
			-b.m_e3_e1, // e3
			b.m_e2_e3 // e1_e2_e3
		);

}
/// <summary>Returns geometric product of I3_t and rotor.
/// </summary>
public static oddVersor gp(I3_t a, rotor b)
{
	return new oddVersor(oddVersor.coord_e1_e2_e3_e1e2e3,
			-b.m_e2_e3, // e1
			-b.m_e3_e1, // e2
			-b.m_e1_e2, // e3
			b.m_scalar // e1_e2_e3
		);

}
/// <summary>Returns geometric product of bivector and e1_t.
/// </summary>
public static oddVersor gp(bivector a, e1_t b)
{
	return new oddVersor(oddVersor.coord_e1_e2_e3_e1e2e3,
			0.0, // e1
			-a.m_e1_e2, // e2
			a.m_e3_e1, // e3
			a.m_e2_e3 // e1_e2_e3
		);

}
/// <summary>Returns a bitmap where each one-bit means that at least one coordinate of the respective grade of  a is larger than b
/// </summary>
public static int gradeBitmap(mv_if a, double b)
{
	int bitmap = 0;
	double[][] ac = a.to_mv().c();
	
	if (ac[0] != null) {
		if (!zeroGroup_0(ac[0], b)) bitmap |= 1;
	}
	
	if (ac[1] != null) {
		if (!zeroGroup_1(ac[1], b)) bitmap |= 2;
	}
	
	if (ac[2] != null) {
		if (!zeroGroup_2(ac[2], b)) bitmap |= 4;
	}
	
	if (ac[3] != null) {
		if (!zeroGroup_3(ac[3], b)) bitmap |= 8;
	}
	return bitmap;
}
/// <summary>Returns a bitmap where each one-bit means that at least one coordinate of the respective grade of  a is larger than b
/// </summary>
public static int gradeBitmap(rotor a, double b)
{
	int bitmap = 0;
	if ((a.m_scalar < -b) || (a.m_scalar > b)) bitmap |= 1;
	if ((a.m_e1_e2 < -b) || (a.m_e1_e2 > b)) bitmap |= 4;
	if ((a.m_e2_e3 < -b) || (a.m_e2_e3 > b)) bitmap |= 4;
	if ((a.m_e3_e1 < -b) || (a.m_e3_e1 > b)) bitmap |= 4;
	return bitmap;
}
/// <summary>Returns a bitmap where each one-bit means that at least one coordinate of the respective grade of  a is larger than b
/// </summary>
public static int gradeBitmap(vector a, double b)
{
	int bitmap = 0;
	if ((a.m_e1 < -b) || (a.m_e1 > b)) bitmap |= 2;
	if ((a.m_e2 < -b) || (a.m_e2 > b)) bitmap |= 2;
	if ((a.m_e3 < -b) || (a.m_e3 > b)) bitmap |= 2;
	return bitmap;
}
/// <summary>Returns a bitmap where each one-bit means that at least one coordinate of the respective grade of  a is larger than b
/// </summary>
public static int gradeBitmap(bivector a, double b)
{
	int bitmap = 0;
	if ((a.m_e1_e2 < -b) || (a.m_e1_e2 > b)) bitmap |= 4;
	if ((a.m_e2_e3 < -b) || (a.m_e2_e3 > b)) bitmap |= 4;
	if ((a.m_e3_e1 < -b) || (a.m_e3_e1 > b)) bitmap |= 4;
	return bitmap;
}
/// <summary>Returns a bitmap where each one-bit means that at least one coordinate of the respective grade of  a is larger than b
/// </summary>
public static int gradeBitmap(trivector a, double b)
{
	int bitmap = 0;
	if ((a.m_e1_e2_e3 < -b) || (a.m_e1_e2_e3 > b)) bitmap |= 8;
	return bitmap;
}
/// <summary>Returns a bitmap where each one-bit means that at least one coordinate of the respective grade of  a is larger than b
/// </summary>
public static int gradeBitmap(e1_t a, double b)
{
	int bitmap = 0;
	if (1.0 > b) bitmap |= 2;
	return bitmap;
}
/// <summary>Returns a bitmap where each one-bit means that at least one coordinate of the respective grade of  a is larger than b
/// </summary>
public static int gradeBitmap(e2_t a, double b)
{
	int bitmap = 0;
	if (1.0 > b) bitmap |= 2;
	return bitmap;
}
/// <summary>Returns a bitmap where each one-bit means that at least one coordinate of the respective grade of  a is larger than b
/// </summary>
public static int gradeBitmap(I3_t a, double b)
{
	int bitmap = 0;
	if (1.0 > b) bitmap |= 8;
	return bitmap;
}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of mv and mv.
/// </summary>
public static mv hp(mv_if a, mv_if b)
{
	double[][] ac = a.to_mv().c();
	double[][] bc = b.to_mv().c();
	double[][] cc = new double[4][];
	
	if (ac[0] != null) {
		if (bc[0] != null) {
			cc[0] = new double[1];
			hp_0_0(ac[0], bc[0], cc[0]);
		}
	}
	
	if (ac[1] != null) {
		if (bc[1] != null) {
			cc[1] = new double[3];
			hp_1_1(ac[1], bc[1], cc[1]);
		}
	}
	
	if (ac[2] != null) {
		if (bc[2] != null) {
			cc[2] = new double[3];
			hp_2_2(ac[2], bc[2], cc[2]);
		}
	}
	
	if (ac[3] != null) {
		if (bc[3] != null) {
			cc[3] = new double[1];
			hp_3_3(ac[3], bc[3], cc[3]);
		}
	}
	return new mv(cc);
}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of vector and vector.
/// </summary>
public static vector hp(vector a, vector b)
{
	return new vector(vector.coord_e1_e2_e3,
			a.m_e1*b.m_e1, // e1
			a.m_e2*b.m_e2, // e2
			a.m_e3*b.m_e3 // e3
		);

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of bivector and bivector.
/// </summary>
public static bivector hp(bivector a, bivector b)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			a.m_e1_e2*b.m_e1_e2, // e1_e2
			a.m_e2_e3*b.m_e2_e3, // e2_e3
			a.m_e3_e1*b.m_e3_e1 // e3_e1
		);

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of rotor and rotor.
/// </summary>
public static rotor hp(rotor a, rotor b)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			a.m_scalar*b.m_scalar, // scalar
			a.m_e1_e2*b.m_e1_e2, // e1_e2
			a.m_e2_e3*b.m_e2_e3, // e2_e3
			a.m_e3_e1*b.m_e3_e1 // e3_e1
		);

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of bivector and rotor.
/// </summary>
public static bivector hp(bivector a, rotor b)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			a.m_e1_e2*b.m_e1_e2, // e1_e2
			a.m_e2_e3*b.m_e2_e3, // e2_e3
			a.m_e3_e1*b.m_e3_e1 // e3_e1
		);

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of trivector and trivector.
/// </summary>
public static trivector hp(trivector a, trivector b)
{
	return new trivector(trivector.coord_e1e2e3,
			a.m_e1_e2_e3*b.m_e1_e2_e3 // e1_e2_e3
		);

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of trivector and oddVersor.
/// </summary>
public static trivector hp(trivector a, oddVersor b)
{
	return new trivector(trivector.coord_e1e2e3,
			a.m_e1_e2_e3*b.m_e1_e2_e3 // e1_e2_e3
		);

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of rotor and bivector.
/// </summary>
public static bivector hp(rotor a, bivector b)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			a.m_e1_e2*b.m_e1_e2, // e1_e2
			a.m_e2_e3*b.m_e2_e3, // e2_e3
			a.m_e3_e1*b.m_e3_e1 // e3_e1
		);

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of e1_t and e1_t.
/// </summary>
public static e1_t hp(e1_t a, e1_t b)
{
	return new e1_t(		);

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of e2_t and e3_t.
/// </summary>
public static double hp(e2_t a, e3_t b)
{
	return 0.0;

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of oddVersor and I3_t.
/// </summary>
public static trivector hp(oddVersor a, I3_t b)
{
	return new trivector(trivector.coord_e1e2e3,
			a.m_e1_e2_e3 // e1_e2_e3
		);

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of mv and mv.
/// </summary>
public static mv ihp(mv_if a, mv_if b)
{
	double[][] ac = a.to_mv().c();
	double[][] bc = b.to_mv().c();
	double[][] cc = new double[4][];
	
	if (ac[0] != null) {
		if (bc[0] != null) {
			cc[0] = new double[1];
			ihp_0_0(ac[0], bc[0], cc[0]);
		}
	}
	
	if (ac[1] != null) {
		if (bc[1] != null) {
			cc[1] = new double[3];
			ihp_1_1(ac[1], bc[1], cc[1]);
		}
	}
	
	if (ac[2] != null) {
		if (bc[2] != null) {
			cc[2] = new double[3];
			ihp_2_2(ac[2], bc[2], cc[2]);
		}
	}
	
	if (ac[3] != null) {
		if (bc[3] != null) {
			cc[3] = new double[1];
			ihp_3_3(ac[3], bc[3], cc[3]);
		}
	}
	return new mv(cc);
}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of vector and vector.
/// </summary>
public static vector ihp(vector a, vector b)
{
	return new vector(vector.coord_e1_e2_e3,
			a.m_e1/((b.m_e1)), // e1
			a.m_e2/((b.m_e2)), // e2
			a.m_e3/((b.m_e3)) // e3
		);

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of bivector and bivector.
/// </summary>
public static bivector ihp(bivector a, bivector b)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			a.m_e1_e2/((b.m_e1_e2)), // e1_e2
			a.m_e2_e3/((b.m_e2_e3)), // e2_e3
			-a.m_e3_e1/((-b.m_e3_e1)) // e3_e1
		);

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of rotor and rotor.
/// </summary>
public static rotor ihp(rotor a, rotor b)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			a.m_scalar/((b.m_scalar)), // scalar
			a.m_e1_e2/((b.m_e1_e2)), // e1_e2
			a.m_e2_e3/((b.m_e2_e3)), // e2_e3
			-a.m_e3_e1/((-b.m_e3_e1)) // e3_e1
		);

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of bivector and rotor.
/// </summary>
public static bivector ihp(bivector a, rotor b)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			a.m_e1_e2/((b.m_e1_e2)), // e1_e2
			a.m_e2_e3/((b.m_e2_e3)), // e2_e3
			-a.m_e3_e1/((-b.m_e3_e1)) // e3_e1
		);

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of rotor and bivector.
/// </summary>
public static bivector ihp(rotor a, bivector b)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			a.m_e1_e2/((b.m_e1_e2)), // e1_e2
			a.m_e2_e3/((b.m_e2_e3)), // e2_e3
			-a.m_e3_e1/((-b.m_e3_e1)) // e3_e1
		);

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of trivector and oddVersor.
/// </summary>
public static trivector ihp(trivector a, oddVersor b)
{
	return new trivector(trivector.coord_e1e2e3,
			a.m_e1_e2_e3/((b.m_e1_e2_e3)) // e1_e2_e3
		);

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of vector and e1_t.
/// </summary>
public static vector ihp(vector a, e1_t b)
{
	return new vector(vector.coord_e1_e2_e3,
			a.m_e1, // e1
			0.0, // e2
			0.0 // e3
		);

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of e2_t and e3_t.
/// </summary>
public static double ihp(e2_t a, e3_t b)
{
	return 0.0;

}
/// <summary>Returns Hadamard product (coordinate-wise multiplication) of trivector and I3_t.
/// </summary>
public static trivector ihp(trivector a, I3_t b)
{
	return new trivector(trivector.coord_e1e2e3,
			a.m_e1_e2_e3 // e1_e2_e3
		);

}
/// <summary>Returns (a + 1).
/// </summary>
public static mv increment(mv_if a)
{
	mv _dst = new mv(a.to_mv());
	double val = _dst.get_scalar() + 1.0;
	_dst.set_scalar(val);
	return _dst;
}
/// <summary>Returns (a + 1).
/// </summary>
public static rotor increment(bivector a)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			1.0, // scalar
			a.m_e1_e2, // e1_e2
			a.m_e2_e3, // e2_e3
			a.m_e3_e1 // e3_e1
		);
}
/// <summary>Returns (a + 1).
/// </summary>
public static rotor increment(rotor a)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			(1.0+a.m_scalar), // scalar
			a.m_e1_e2, // e1_e2
			a.m_e2_e3, // e2_e3
			a.m_e3_e1 // e3_e1
		);
}
/// <summary>Returns (a - 1).
/// </summary>
public static mv decrement(mv_if a)
{
	mv _dst = new mv(a.to_mv());
	double val = _dst.get_scalar() - 1.0;
	_dst.set_scalar(val);
	return _dst;
}
/// <summary>Returns (a - 1).
/// </summary>
public static rotor decrement(bivector a)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			-1.0, // scalar
			a.m_e1_e2, // e1_e2
			a.m_e2_e3, // e2_e3
			a.m_e3_e1 // e3_e1
		);
}
/// <summary>Returns (a - 1).
/// </summary>
public static rotor decrement(rotor a)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			(-1.0+a.m_scalar), // scalar
			a.m_e1_e2, // e1_e2
			a.m_e2_e3, // e2_e3
			a.m_e3_e1 // e3_e1
		);
}
/// <summary>Returns scalar product of mv and mv.
/// </summary>
public static double sp(mv_if a, mv_if b)
{
	double[][] ac = a.to_mv().c();
	double[][] bc = b.to_mv().c();
	double[][] cc = new double[4][];
	cc[0] = new double[1];
	if (ac[0] != null) {
		if (bc[0] != null) {
			gp_default_0_0_0(ac[0], bc[0], cc[0]);
		}
	}
	if (ac[1] != null) {
		if (bc[1] != null) {
			gp_default_1_1_0(ac[1], bc[1], cc[0]);
		}
	}
	if (ac[2] != null) {
		if (bc[2] != null) {
			gp_default_2_2_0(ac[2], bc[2], cc[0]);
		}
	}
	if (ac[3] != null) {
		if (bc[3] != null) {
			gp_default_3_3_0(ac[3], bc[3], cc[0]);
		}
	}
	return cc[0][0];
}
/// <summary>Returns left contraction of mv and mv.
/// </summary>
public static mv lc(mv_if a, mv_if b)
{
	double[][] ac = a.to_mv().c();
	double[][] bc = b.to_mv().c();
	double[][] cc = new double[4][];
	if (ac[0] != null) {
		if (bc[0] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_0_0_0(ac[0], bc[0], cc[0]);
		}
		if (bc[1] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_0_1_1(ac[0], bc[1], cc[1]);
		}
		if (bc[2] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_0_2_2(ac[0], bc[2], cc[2]);
		}
		if (bc[3] != null) {
			if (cc[3] == null) cc[3] = new double[1];
			gp_default_0_3_3(ac[0], bc[3], cc[3]);
		}
	}
	if (ac[1] != null) {
		if (bc[1] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_1_1_0(ac[1], bc[1], cc[0]);
		}
		if (bc[2] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_1_2_1(ac[1], bc[2], cc[1]);
		}
		if (bc[3] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_1_3_2(ac[1], bc[3], cc[2]);
		}
	}
	if (ac[2] != null) {
		if (bc[2] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_2_2_0(ac[2], bc[2], cc[0]);
		}
		if (bc[3] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_2_3_1(ac[2], bc[3], cc[1]);
		}
	}
	if (ac[3] != null) {
		if (bc[3] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_3_3_0(ac[3], bc[3], cc[0]);
		}
	}
	return new mv(cc);
}
/// <summary>Returns right contraction of mv and mv.
/// </summary>
public static mv rc(mv_if a, mv_if b)
{
	double[][] ac = a.to_mv().c();
	double[][] bc = b.to_mv().c();
	double[][] cc = new double[4][];
	if (ac[0] != null) {
		if (bc[0] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_0_0_0(ac[0], bc[0], cc[0]);
		}
	}
	if (ac[1] != null) {
		if (bc[0] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_1_0_1(ac[1], bc[0], cc[1]);
		}
		if (bc[1] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_1_1_0(ac[1], bc[1], cc[0]);
		}
	}
	if (ac[2] != null) {
		if (bc[0] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_2_0_2(ac[2], bc[0], cc[2]);
		}
		if (bc[1] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_2_1_1(ac[2], bc[1], cc[1]);
		}
		if (bc[2] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_2_2_0(ac[2], bc[2], cc[0]);
		}
	}
	if (ac[3] != null) {
		if (bc[0] != null) {
			if (cc[3] == null) cc[3] = new double[1];
			gp_default_3_0_3(ac[3], bc[0], cc[3]);
		}
		if (bc[1] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_3_1_2(ac[3], bc[1], cc[2]);
		}
		if (bc[2] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_3_2_1(ac[3], bc[2], cc[1]);
		}
		if (bc[3] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_3_3_0(ac[3], bc[3], cc[0]);
		}
	}
	return new mv(cc);
}
/// <summary>Returns Hestenes inner product of mv and mv.
/// </summary>
public static mv hip(mv_if a, mv_if b)
{
	double[][] ac = a.to_mv().c();
	double[][] bc = b.to_mv().c();
	double[][] cc = new double[4][];
	if (ac[1] != null) {
		if (bc[1] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_1_1_0(ac[1], bc[1], cc[0]);
		}
		if (bc[2] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_1_2_1(ac[1], bc[2], cc[1]);
		}
		if (bc[3] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_1_3_2(ac[1], bc[3], cc[2]);
		}
	}
	if (ac[2] != null) {
		if (bc[1] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_2_1_1(ac[2], bc[1], cc[1]);
		}
		if (bc[2] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_2_2_0(ac[2], bc[2], cc[0]);
		}
		if (bc[3] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_2_3_1(ac[2], bc[3], cc[1]);
		}
	}
	if (ac[3] != null) {
		if (bc[1] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_3_1_2(ac[3], bc[1], cc[2]);
		}
		if (bc[2] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_3_2_1(ac[3], bc[2], cc[1]);
		}
		if (bc[3] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_3_3_0(ac[3], bc[3], cc[0]);
		}
	}
	return new mv(cc);
}
/// <summary>Returns Modified Hestenes inner product of mv and mv.
/// </summary>
public static mv mhip(mv_if a, mv_if b)
{
	double[][] ac = a.to_mv().c();
	double[][] bc = b.to_mv().c();
	double[][] cc = new double[4][];
	if (ac[0] != null) {
		if (bc[0] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_0_0_0(ac[0], bc[0], cc[0]);
		}
		if (bc[1] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_0_1_1(ac[0], bc[1], cc[1]);
		}
		if (bc[2] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_0_2_2(ac[0], bc[2], cc[2]);
		}
		if (bc[3] != null) {
			if (cc[3] == null) cc[3] = new double[1];
			gp_default_0_3_3(ac[0], bc[3], cc[3]);
		}
	}
	if (ac[1] != null) {
		if (bc[0] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_1_0_1(ac[1], bc[0], cc[1]);
		}
		if (bc[1] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_1_1_0(ac[1], bc[1], cc[0]);
		}
		if (bc[2] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_1_2_1(ac[1], bc[2], cc[1]);
		}
		if (bc[3] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_1_3_2(ac[1], bc[3], cc[2]);
		}
	}
	if (ac[2] != null) {
		if (bc[0] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_2_0_2(ac[2], bc[0], cc[2]);
		}
		if (bc[1] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_2_1_1(ac[2], bc[1], cc[1]);
		}
		if (bc[2] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_2_2_0(ac[2], bc[2], cc[0]);
		}
		if (bc[3] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_2_3_1(ac[2], bc[3], cc[1]);
		}
	}
	if (ac[3] != null) {
		if (bc[0] != null) {
			if (cc[3] == null) cc[3] = new double[1];
			gp_default_3_0_3(ac[3], bc[0], cc[3]);
		}
		if (bc[1] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_3_1_2(ac[3], bc[1], cc[2]);
		}
		if (bc[2] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_3_2_1(ac[3], bc[2], cc[1]);
		}
		if (bc[3] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_3_3_0(ac[3], bc[3], cc[0]);
		}
	}
	return new mv(cc);
}
/// <summary>Returns scalar product of vector and vector.
/// </summary>
public static double sp(vector a, vector b)
{
	return (a.m_e1*b.m_e1+a.m_e2*b.m_e2+a.m_e3*b.m_e3);

}
/// <summary>Returns left contraction of vector and vector.
/// </summary>
public static double lc(vector a, vector b)
{
	return (a.m_e1*b.m_e1+a.m_e2*b.m_e2+a.m_e3*b.m_e3);

}
/// <summary>Returns right contraction of vector and vector.
/// </summary>
public static double rc(vector a, vector b)
{
	return (a.m_e1*b.m_e1+a.m_e2*b.m_e2+a.m_e3*b.m_e3);

}
/// <summary>Returns Hestenes inner product of vector and vector.
/// </summary>
public static double hip(vector a, vector b)
{
	return (a.m_e1*b.m_e1+a.m_e2*b.m_e2+a.m_e3*b.m_e3);

}
/// <summary>Returns Modified Hestenes inner product of vector and vector.
/// </summary>
public static double mhip(vector a, vector b)
{
	return (a.m_e1*b.m_e1+a.m_e2*b.m_e2+a.m_e3*b.m_e3);

}
/// <summary>Returns scalar product of bivector and vector.
/// </summary>
public static double sp(bivector a, vector b)
{
	return 0.0;

}
/// <summary>Returns left contraction of bivector and vector.
/// </summary>
public static double lc(bivector a, vector b)
{
	return 0.0;

}
/// <summary>Returns right contraction of bivector and vector.
/// </summary>
public static vector rc(bivector a, vector b)
{
	return new vector(vector.coord_e1_e2_e3,
			(a.m_e1_e2*b.m_e2-a.m_e3_e1*b.m_e3), // e1
			(-a.m_e1_e2*b.m_e1+a.m_e2_e3*b.m_e3), // e2
			(-a.m_e2_e3*b.m_e2+a.m_e3_e1*b.m_e1) // e3
		);

}
/// <summary>Returns Hestenes inner product of bivector and vector.
/// </summary>
public static vector hip(bivector a, vector b)
{
	return new vector(vector.coord_e1_e2_e3,
			(a.m_e1_e2*b.m_e2-a.m_e3_e1*b.m_e3), // e1
			(-a.m_e1_e2*b.m_e1+a.m_e2_e3*b.m_e3), // e2
			(-a.m_e2_e3*b.m_e2+a.m_e3_e1*b.m_e1) // e3
		);

}
/// <summary>Returns Modified Hestenes inner product of bivector and vector.
/// </summary>
public static vector mhip(bivector a, vector b)
{
	return new vector(vector.coord_e1_e2_e3,
			(a.m_e1_e2*b.m_e2-a.m_e3_e1*b.m_e3), // e1
			(-a.m_e1_e2*b.m_e1+a.m_e2_e3*b.m_e3), // e2
			(-a.m_e2_e3*b.m_e2+a.m_e3_e1*b.m_e1) // e3
		);

}
/// <summary>Returns scalar product of trivector and trivector.
/// </summary>
public static double sp(trivector a, trivector b)
{
	return -a.m_e1_e2_e3*b.m_e1_e2_e3;

}
/// <summary>Returns left contraction of trivector and trivector.
/// </summary>
public static double lc(trivector a, trivector b)
{
	return -a.m_e1_e2_e3*b.m_e1_e2_e3;

}
/// <summary>Returns right contraction of trivector and trivector.
/// </summary>
public static double rc(trivector a, trivector b)
{
	return -a.m_e1_e2_e3*b.m_e1_e2_e3;

}
/// <summary>Returns Hestenes inner product of trivector and trivector.
/// </summary>
public static double hip(trivector a, trivector b)
{
	return -a.m_e1_e2_e3*b.m_e1_e2_e3;

}
/// <summary>Returns Modified Hestenes inner product of trivector and trivector.
/// </summary>
public static double mhip(trivector a, trivector b)
{
	return -a.m_e1_e2_e3*b.m_e1_e2_e3;

}
/// <summary>Returns scalar product of vector and bivector.
/// </summary>
public static double sp(vector a, bivector b)
{
	return 0.0;

}
/// <summary>Returns left contraction of vector and bivector.
/// </summary>
public static vector lc(vector a, bivector b)
{
	return new vector(vector.coord_e1_e2_e3,
			(-a.m_e2*b.m_e1_e2+a.m_e3*b.m_e3_e1), // e1
			(a.m_e1*b.m_e1_e2-a.m_e3*b.m_e2_e3), // e2
			(-a.m_e1*b.m_e3_e1+a.m_e2*b.m_e2_e3) // e3
		);

}
/// <summary>Returns right contraction of vector and bivector.
/// </summary>
public static double rc(vector a, bivector b)
{
	return 0.0;

}
/// <summary>Returns Hestenes inner product of vector and bivector.
/// </summary>
public static vector hip(vector a, bivector b)
{
	return new vector(vector.coord_e1_e2_e3,
			(-a.m_e2*b.m_e1_e2+a.m_e3*b.m_e3_e1), // e1
			(a.m_e1*b.m_e1_e2-a.m_e3*b.m_e2_e3), // e2
			(-a.m_e1*b.m_e3_e1+a.m_e2*b.m_e2_e3) // e3
		);

}
/// <summary>Returns Modified Hestenes inner product of vector and bivector.
/// </summary>
public static vector mhip(vector a, bivector b)
{
	return new vector(vector.coord_e1_e2_e3,
			(-a.m_e2*b.m_e1_e2+a.m_e3*b.m_e3_e1), // e1
			(a.m_e1*b.m_e1_e2-a.m_e3*b.m_e2_e3), // e2
			(-a.m_e1*b.m_e3_e1+a.m_e2*b.m_e2_e3) // e3
		);

}
/// <summary>Returns scalar product of vector and rotor.
/// </summary>
public static double sp(vector a, rotor b)
{
	return 0.0;

}
/// <summary>Returns left contraction of vector and rotor.
/// </summary>
public static vector lc(vector a, rotor b)
{
	return new vector(vector.coord_e1_e2_e3,
			(-a.m_e2*b.m_e1_e2+a.m_e3*b.m_e3_e1), // e1
			(a.m_e1*b.m_e1_e2-a.m_e3*b.m_e2_e3), // e2
			(-a.m_e1*b.m_e3_e1+a.m_e2*b.m_e2_e3) // e3
		);

}
/// <summary>Returns right contraction of vector and rotor.
/// </summary>
public static vector rc(vector a, rotor b)
{
	return new vector(vector.coord_e1_e2_e3,
			a.m_e1*b.m_scalar, // e1
			a.m_e2*b.m_scalar, // e2
			a.m_e3*b.m_scalar // e3
		);

}
/// <summary>Returns Hestenes inner product of vector and rotor.
/// </summary>
public static vector hip(vector a, rotor b)
{
	return new vector(vector.coord_e1_e2_e3,
			(-a.m_e2*b.m_e1_e2+a.m_e3*b.m_e3_e1), // e1
			(a.m_e1*b.m_e1_e2-a.m_e3*b.m_e2_e3), // e2
			(-a.m_e1*b.m_e3_e1+a.m_e2*b.m_e2_e3) // e3
		);

}
/// <summary>Returns Modified Hestenes inner product of vector and rotor.
/// </summary>
public static vector mhip(vector a, rotor b)
{
	return new vector(vector.coord_e1_e2_e3,
			(a.m_e1*b.m_scalar-a.m_e2*b.m_e1_e2+a.m_e3*b.m_e3_e1), // e1
			(a.m_e1*b.m_e1_e2+a.m_e2*b.m_scalar-a.m_e3*b.m_e2_e3), // e2
			(-a.m_e1*b.m_e3_e1+a.m_e2*b.m_e2_e3+a.m_e3*b.m_scalar) // e3
		);

}
/// <summary>Returns scalar product of rotor and bivector.
/// </summary>
public static double sp(rotor a, bivector b)
{
	return (-a.m_e1_e2*b.m_e1_e2-a.m_e2_e3*b.m_e2_e3-a.m_e3_e1*b.m_e3_e1);

}
/// <summary>Returns left contraction of rotor and bivector.
/// </summary>
public static rotor lc(rotor a, bivector b)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			(-a.m_e1_e2*b.m_e1_e2-a.m_e2_e3*b.m_e2_e3-a.m_e3_e1*b.m_e3_e1), // scalar
			a.m_scalar*b.m_e1_e2, // e1_e2
			a.m_scalar*b.m_e2_e3, // e2_e3
			a.m_scalar*b.m_e3_e1 // e3_e1
		);

}
/// <summary>Returns right contraction of rotor and bivector.
/// </summary>
public static double rc(rotor a, bivector b)
{
	return (-a.m_e1_e2*b.m_e1_e2-a.m_e2_e3*b.m_e2_e3-a.m_e3_e1*b.m_e3_e1);

}
/// <summary>Returns Hestenes inner product of rotor and bivector.
/// </summary>
public static double hip(rotor a, bivector b)
{
	return (-a.m_e1_e2*b.m_e1_e2-a.m_e2_e3*b.m_e2_e3-a.m_e3_e1*b.m_e3_e1);

}
/// <summary>Returns Modified Hestenes inner product of rotor and bivector.
/// </summary>
public static rotor mhip(rotor a, bivector b)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			(-a.m_e1_e2*b.m_e1_e2-a.m_e2_e3*b.m_e2_e3-a.m_e3_e1*b.m_e3_e1), // scalar
			a.m_scalar*b.m_e1_e2, // e1_e2
			a.m_scalar*b.m_e2_e3, // e2_e3
			a.m_scalar*b.m_e3_e1 // e3_e1
		);

}
/// <summary>Returns scalar product of rotor and trivector.
/// </summary>
public static double sp(rotor a, trivector b)
{
	return 0.0;

}
/// <summary>Returns left contraction of rotor and trivector.
/// </summary>
public static oddVersor lc(rotor a, trivector b)
{
	return new oddVersor(oddVersor.coord_e1_e2_e3_e1e2e3,
			-a.m_e2_e3*b.m_e1_e2_e3, // e1
			-a.m_e3_e1*b.m_e1_e2_e3, // e2
			-a.m_e1_e2*b.m_e1_e2_e3, // e3
			a.m_scalar*b.m_e1_e2_e3 // e1_e2_e3
		);

}
/// <summary>Returns right contraction of rotor and trivector.
/// </summary>
public static double rc(rotor a, trivector b)
{
	return 0.0;

}
/// <summary>Returns Hestenes inner product of rotor and trivector.
/// </summary>
public static vector hip(rotor a, trivector b)
{
	return new vector(vector.coord_e1_e2_e3,
			-a.m_e2_e3*b.m_e1_e2_e3, // e1
			-a.m_e3_e1*b.m_e1_e2_e3, // e2
			-a.m_e1_e2*b.m_e1_e2_e3 // e3
		);

}
/// <summary>Returns Modified Hestenes inner product of rotor and trivector.
/// </summary>
public static oddVersor mhip(rotor a, trivector b)
{
	return new oddVersor(oddVersor.coord_e1_e2_e3_e1e2e3,
			-a.m_e2_e3*b.m_e1_e2_e3, // e1
			-a.m_e3_e1*b.m_e1_e2_e3, // e2
			-a.m_e1_e2*b.m_e1_e2_e3, // e3
			a.m_scalar*b.m_e1_e2_e3 // e1_e2_e3
		);

}
/// <summary>Returns scalar product of e1_t and I3_t.
/// </summary>
public static double sp(e1_t a, I3_t b)
{
	return 0.0;

}
/// <summary>Returns left contraction of I3_t and e3_t.
/// </summary>
public static double lc(I3_t a, e3_t b)
{
	return 0.0;

}
/// <summary>Returns right contraction of e1_t and e1_t.
/// </summary>
public static double rc(e1_t a, e1_t b)
{
	return 1.0;

}
/// <summary>Returns Hestenes inner product of e2_t and I3_t.
/// </summary>
public static bivector hip(e2_t a, I3_t b)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			0.0, // e1_e2
			0.0, // e2_e3
			1.0 // e3_e1
		);

}
/// <summary>Returns Modified Hestenes inner product of I3_t and I3_t.
/// </summary>
public static double mhip(I3_t a, I3_t b)
{
	return -1.0;

}
/// <summary>Returns norm of mv using default metric.
/// </summary>
public static double norm(mv_if a)
{
	double[] c = new double[1];
	double[][] ac = a.to_mv().c();
	double n2 = 0.0;
	
	if (ac[0] != null) { /* group 0 (grade 0) */
		c[0] = 0.0;
		gp_default_0_0_0(ac[0], ac[ 0], c);
		n2 += c[0];
	}
	
	if (ac[1] != null) { /* group 1 (grade 1) */
		c[0] = 0.0;
		gp_default_1_1_0(ac[1], ac[ 1], c);
		n2 += c[0];
	}
	
	if (ac[2] != null) { /* group 2 (grade 2) */
		c[0] = 0.0;
		gp_default_2_2_0(ac[2], ac[ 2], c);
		n2 -= c[0];
	}
	
	if (ac[3] != null) { /* group 3 (grade 3) */
		c[0] = 0.0;
		gp_default_3_3_0(ac[3], ac[ 3], c);
		n2 -= c[0];
	}
	return Math.Sqrt(n2);
}
/// <summary>internal conversion function
/// </summary>
public static double norm_returns_scalar(mv a) {
	return norm(a);
}
/// <summary>Returns norm of vector using default metric.
/// </summary>
public static double norm(vector a)
{
	return Math.Abs(Math.Sqrt((a.m_e1*a.m_e1+a.m_e2*a.m_e2+a.m_e3*a.m_e3)));

}
/// <summary>internal conversion function
/// </summary>
public static double norm_returns_scalar(vector a) {
	return norm(a);
}
/// <summary>Returns norm of bivector using default metric.
/// </summary>
public static double norm(bivector a)
{
	return Math.Abs(Math.Sqrt((a.m_e1_e2*a.m_e1_e2+a.m_e2_e3*a.m_e2_e3+a.m_e3_e1*a.m_e3_e1)));

}
/// <summary>internal conversion function
/// </summary>
public static double norm_returns_scalar(bivector a) {
	return norm(a);
}
/// <summary>Returns norm of trivector using default metric.
/// </summary>
public static double norm(trivector a)
{
	return Math.Abs(Math.Sqrt(a.m_e1_e2_e3*a.m_e1_e2_e3));

}
/// <summary>internal conversion function
/// </summary>
public static double norm_returns_scalar(trivector a) {
	return norm(a);
}
/// <summary>Returns norm of rotor using default metric.
/// </summary>
public static double norm(rotor a)
{
	return Math.Abs(Math.Sqrt((a.m_e1_e2*a.m_e1_e2+a.m_e2_e3*a.m_e2_e3+a.m_e3_e1*a.m_e3_e1+a.m_scalar*a.m_scalar)));

}
/// <summary>internal conversion function
/// </summary>
public static double norm_returns_scalar(rotor a) {
	return norm(a);
}
/// <summary>Returns norm of e1_t using default metric.
/// </summary>
public static double norm(e1_t a)
{
	return Math.Abs(1.0);

}
/// <summary>internal conversion function
/// </summary>
public static double norm_returns_scalar(e1_t a) {
	return norm(a);
}
/// <summary>Returns norm of e3_t using default metric.
/// </summary>
public static double norm(e3_t a)
{
	return Math.Abs(1.0);

}
/// <summary>internal conversion function
/// </summary>
public static double norm_returns_scalar(e3_t a) {
	return norm(a);
}
/// <summary>Returns norm of I3_t using default metric.
/// </summary>
public static double norm(I3_t a)
{
	return Math.Abs(1.0);

}
/// <summary>internal conversion function
/// </summary>
public static double norm_returns_scalar(I3_t a) {
	return norm(a);
}
/// <summary>Returns norm2 of mv using default metric.
/// </summary>
public static double norm2(mv_if a)
{
	double[] c = new double[1];
	double[][] ac = a.to_mv().c();
	double n2 = 0.0;
	
	if (ac[0] != null) { /* group 0 (grade 0) */
		c[0] = 0.0;
		gp_default_0_0_0(ac[0], ac[ 0], c);
		n2 += c[0];
	}
	
	if (ac[1] != null) { /* group 1 (grade 1) */
		c[0] = 0.0;
		gp_default_1_1_0(ac[1], ac[ 1], c);
		n2 += c[0];
	}
	
	if (ac[2] != null) { /* group 2 (grade 2) */
		c[0] = 0.0;
		gp_default_2_2_0(ac[2], ac[ 2], c);
		n2 -= c[0];
	}
	
	if (ac[3] != null) { /* group 3 (grade 3) */
		c[0] = 0.0;
		gp_default_3_3_0(ac[3], ac[ 3], c);
		n2 -= c[0];
	}
	return n2;
}
/// <summary>internal conversion function
/// </summary>
public static double norm2_returns_scalar(mv a) {
	return norm2(a);
}
/// <summary>Returns norm2 of vector using default metric.
/// </summary>
public static double norm2(vector a)
{
	return (a.m_e1*a.m_e1+a.m_e2*a.m_e2+a.m_e3*a.m_e3);

}
/// <summary>internal conversion function
/// </summary>
public static double norm2_returns_scalar(vector a) {
	return norm2(a);
}
/// <summary>Returns norm2 of bivector using default metric.
/// </summary>
public static double norm2(bivector a)
{
	return (a.m_e1_e2*a.m_e1_e2+a.m_e2_e3*a.m_e2_e3+a.m_e3_e1*a.m_e3_e1);

}
/// <summary>internal conversion function
/// </summary>
public static double norm2_returns_scalar(bivector a) {
	return norm2(a);
}
/// <summary>Returns norm2 of trivector using default metric.
/// </summary>
public static double norm2(trivector a)
{
	return a.m_e1_e2_e3*a.m_e1_e2_e3;

}
/// <summary>internal conversion function
/// </summary>
public static double norm2_returns_scalar(trivector a) {
	return norm2(a);
}
/// <summary>Returns norm2 of rotor using default metric.
/// </summary>
public static double norm2(rotor a)
{
	return (a.m_e1_e2*a.m_e1_e2+a.m_e2_e3*a.m_e2_e3+a.m_e3_e1*a.m_e3_e1+a.m_scalar*a.m_scalar);

}
/// <summary>internal conversion function
/// </summary>
public static double norm2_returns_scalar(rotor a) {
	return norm2(a);
}
/// <summary>Returns norm2 of e1_t using default metric.
/// </summary>
public static double norm2(e1_t a)
{
	return 1.0;

}
/// <summary>internal conversion function
/// </summary>
public static double norm2_returns_scalar(e1_t a) {
	return norm2(a);
}
/// <summary>Returns norm2 of e3_t using default metric.
/// </summary>
public static double norm2(e3_t a)
{
	return 1.0;

}
/// <summary>internal conversion function
/// </summary>
public static double norm2_returns_scalar(e3_t a) {
	return norm2(a);
}
/// <summary>Returns norm2 of I3_t using default metric.
/// </summary>
public static double norm2(I3_t a)
{
	return 1.0;

}
/// <summary>internal conversion function
/// </summary>
public static double norm2_returns_scalar(I3_t a) {
	return norm2(a);
}
/// <summary>Returns outer product of mv and mv.
/// </summary>
public static mv op(mv_if a, mv_if b)
{
	double[][] ac = a.to_mv().c();
	double[][] bc = b.to_mv().c();
	double[][] cc = new double[4][];
	if (ac[0] != null) {
		if (bc[0] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_0_0_0(ac[0], bc[0], cc[0]);
		}
		if (bc[1] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_0_1_1(ac[0], bc[1], cc[1]);
		}
		if (bc[2] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_0_2_2(ac[0], bc[2], cc[2]);
		}
		if (bc[3] != null) {
			if (cc[3] == null) cc[3] = new double[1];
			gp_default_0_3_3(ac[0], bc[3], cc[3]);
		}
	}
	if (ac[1] != null) {
		if (bc[0] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_1_0_1(ac[1], bc[0], cc[1]);
		}
		if (bc[1] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_1_1_2(ac[1], bc[1], cc[2]);
		}
		if (bc[2] != null) {
			if (cc[3] == null) cc[3] = new double[1];
			gp_default_1_2_3(ac[1], bc[2], cc[3]);
		}
	}
	if (ac[2] != null) {
		if (bc[0] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_2_0_2(ac[2], bc[0], cc[2]);
		}
		if (bc[1] != null) {
			if (cc[3] == null) cc[3] = new double[1];
			gp_default_2_1_3(ac[2], bc[1], cc[3]);
		}
	}
	if (ac[3] != null) {
		if (bc[0] != null) {
			if (cc[3] == null) cc[3] = new double[1];
			gp_default_3_0_3(ac[3], bc[0], cc[3]);
		}
	}
	return new mv(cc);
}
/// <summary>Returns outer product of vector and vector.
/// </summary>
public static bivector op(vector a, vector b)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			(a.m_e1*b.m_e2-a.m_e2*b.m_e1), // e1_e2
			(a.m_e2*b.m_e3-a.m_e3*b.m_e2), // e2_e3
			-(a.m_e1*b.m_e3-a.m_e3*b.m_e1) // e3_e1
		);

}
/// <summary>Returns outer product of bivector and bivector.
/// </summary>
public static double op(bivector a, bivector b)
{
	return 0.0;

}
/// <summary>Returns outer product of vector and rotor.
/// </summary>
public static oddVersor op(vector a, rotor b)
{
	return new oddVersor(oddVersor.coord_e1_e2_e3_e1e2e3,
			a.m_e1*b.m_scalar, // e1
			a.m_e2*b.m_scalar, // e2
			a.m_e3*b.m_scalar, // e3
			(a.m_e1*b.m_e2_e3+a.m_e2*b.m_e3_e1+a.m_e3*b.m_e1_e2) // e1_e2_e3
		);

}
/// <summary>Returns outer product of vector and trivector.
/// </summary>
public static double op(vector a, trivector b)
{
	return 0.0;

}
/// <summary>Returns outer product of e1_t and e1_t.
/// </summary>
public static double op(e1_t a, e1_t b)
{
	return 0.0;

}
/// <summary>Returns outer product of e1_t and e2_t.
/// </summary>
public static bivector op(e1_t a, e2_t b)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			1.0, // e1_e2
			0.0, // e2_e3
			0.0 // e3_e1
		);

}
/// <summary>Returns outer product of e2_t and e3_t.
/// </summary>
public static bivector op(e2_t a, e3_t b)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			0.0, // e1_e2
			1.0, // e2_e3
			0.0 // e3_e1
		);

}
/// <summary>Returns geometric product of mv and double.
/// </summary>
public static mv gp(mv_if a, double b)
{
	double[][] ac = a.to_mv().c();
	double[][] bc = new double[][]{new double[]{b}};
	double[][] cc = new double[4][];
	if (ac[0] != null) {
			if (cc[0] == null) cc[0] = new double[1];
			gp_default_0_0_0(ac[0], bc[0], cc[0]);
	}
	if (ac[1] != null) {
			if (cc[1] == null) cc[1] = new double[3];
			gp_default_1_0_1(ac[1], bc[0], cc[1]);
	}
	if (ac[2] != null) {
			if (cc[2] == null) cc[2] = new double[3];
			gp_default_2_0_2(ac[2], bc[0], cc[2]);
	}
	if (ac[3] != null) {
			if (cc[3] == null) cc[3] = new double[1];
			gp_default_3_0_3(ac[3], bc[0], cc[3]);
	}
	return new mv(cc);
}
/// <summary>Returns double b * mv a + double c.
/// </summary>
public static mv sas(mv_if a, double b, double c)
{
	double[][] ac = a.to_mv().c();
	double[][] cc = new double[4][];
	
	if (ac[0] != null) {
		cc[0] = new double[1];
		copyMul_0(ac[0], cc[0], b);
		cc[0][0] += c;
	}
	else if (c != 0.0) {
		cc[0] = new double[1];
	cc[0][0] = c;
	}
	
	if (ac[1] != null) {
		cc[1] = new double[3];
		copyMul_1(ac[1], cc[1], b);
	}
	
	if (ac[2] != null) {
		cc[2] = new double[3];
		copyMul_2(ac[2], cc[2], b);
	}
	
	if (ac[3] != null) {
		cc[3] = new double[1];
		copyMul_3(ac[3], cc[3], b);
	}
	return new mv(cc);
}

/// <summary>Computes exponential of mv up to 12th term.
/// 
/// </summary>
public static mv exp(mv x) {
	return exp(x, 12);
}

/// <summary>Computes exponential of mv.
/// 
/// </summary>
public static mv exp(mv x, int order) {
   
	{ // First try special cases: check if (x * x) is scalar
		mv xSquared = gp(x, x);
		double s_xSquared = xSquared.get_scalar();
		if ((norm2_returns_scalar(xSquared) - s_xSquared * s_xSquared) < 1E-14) {
			// OK (x * x == ~scalar), so use special cases:
			if (s_xSquared < 0.0) {
				double a = Math.Sqrt(-s_xSquared);
				return sas(x, Math.Sin(a) / a, Math.Cos(a));
			}
			else if (s_xSquared > 0.0) {
				double a = Math.Sqrt(s_xSquared);
				return sas(x, Math.Sinh(a) / a, Math.Cosh(a));
			}
			else {
				return sas(x, 1.0, 1.0);
			}
		}
	}

	// else do general series eval . . .

	// result = 1 + ....	
	mv result = new mv(1.0);
	if (order == 0) return result;

	// find scale (power of 2) such that its norm is < 1
	ulong maxC = (ulong)x.LargestCoordinate();
	int scale = 1;
	if (maxC > 1) scale <<= 1;
	while (maxC != 0)
	{
		maxC >>= 1;
		scale <<= 1;
	}

	// scale
	mv xScaled = gp(x, 1.0 / (double)scale); 

	// taylor series approximation
	mv xPow1 = new mv(1.0); 
	for (int i = 1; i <= order; i++) {
		mv xPow2 = gp(xPow1, xScaled);
		xPow1 = gp(xPow2, 1.0 / (double)i);
		
		result = add(result, xPow1); // result2 = result1 + xPow1
    }

	// undo scaling
	while (scale > 1)
	{
		result = gp(result, result);
		scale >>= 1;
	}
    
    return result;
} // end of exp()

/// <summary>exp of bivector (uses fast special case)
/// </summary>
public static rotor exp(bivector a)
{
	double _alpha = Math.Sqrt(Math.Abs((-a.m_e1_e2*a.m_e1_e2-a.m_e2_e3*a.m_e2_e3-a.m_e3_e1*a.m_e3_e1)));

	double _mul;
	if (_alpha != 0.0) {
		_mul = Math.Sin(_alpha)/((_alpha));

	}
	else {
		_mul = 0.0;

	}
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			Math.Cos(_alpha), // scalar
			_mul*a.m_e1_e2, // e1_e2
			_mul*a.m_e2_e3, // e2_e3
			_mul*a.m_e3_e1 // e3_e1
		);
}

/// <summary>Computes hyperbolic cosine of mv up to 12th term.
/// 
/// </summary>
public static mv cosh(mv x) {
	return cosh(x, 12);
}

/// <summary>Computes hyperbolic cosine of mv.
/// 
/// </summary>
public static mv cosh(mv x, int order) {
   
	{ // First try special cases: check if (x * x) is scalar
		mv xSquared = gp(x, x);
		double s_xSquared = xSquared.get_scalar();
		if ((norm2_returns_scalar(xSquared) - s_xSquared * s_xSquared) < 1E-14) {
			// OK (x * x == ~scalar), so use special cases:
			if (s_xSquared > 0.0) {
				return new mv(Math.Cosh(Math.Sqrt(s_xSquared)));
			}
			else if (s_xSquared < 0.0) {
				return new mv(Math.Cos(Math.Sqrt(-s_xSquared)));
			}
			else {
				return new mv(1.0);
			}
		}
	}

	// else do general series eval . . .


	mv result = new mv(1.0);
	if (order == 0) return result;

	// taylor series approximation
	mv xPow1 = new mv(1.0);
	for (int i = 1; i <= order; i++) {
		mv xPow2 = gp(xPow1, x);
		xPow1 = gp(xPow2, 1.0 / (double)i); // xPow1 = xScaled^i / i! 
		
		if ((i % 2) == 0) {
			result = add(result, xPow1); 
		}
    }

    return result;
} // end of cosh()
/// <summary>cosh of bivector (uses fast special case)
/// </summary>
public static double cosh(bivector a)
{
	double _alpha = Math.Sqrt(Math.Abs((-a.m_e1_e2*a.m_e1_e2-a.m_e2_e3*a.m_e2_e3-a.m_e3_e1*a.m_e3_e1)));

	return Math.Cos(_alpha);
}

/// <summary>Computes hyperbolic sine of mv up to 12th term.
/// 
/// </summary>
public static mv sinh(mv x) {
	return sinh(x, 12);
}

/// <summary>Computes hyperbolic sine of mv.
/// 
/// </summary>
public static mv sinh(mv x, int order) {
   
	{ // First try special cases: check if (x * x) is scalar
		mv xSquared = gp(x, x);
		double s_xSquared = xSquared.get_scalar();
		if ((norm2_returns_scalar(xSquared) - s_xSquared * s_xSquared) < 1E-14) {
			// OK (x * x == ~scalar), so use special cases:
			if (s_xSquared < 0.0) {
				double a = Math.Sqrt(-s_xSquared);
				return sas(x, Math.Sin(a) / a, 0.0);
			}
			else if (s_xSquared > 0.0) {
				double a = Math.Sqrt(s_xSquared);
				return sas(x, Math.Sinh(a) / a, 0.0);
			}
			else {
				return x;
			}
		}
	}

	// else do general series eval . . .

	// result = A +  A^3/3! + A^5/5!
	mv result = new mv(); // result = 0
    if (order == 0) return result;
    	
	// taylor series approximation
	mv xPow1 = new mv(1.0);
	for (int i = 1; i <= order; i++) {
		mv xPow2 = gp(xPow1, x); // xPow2 = xPow1 * x
		xPow1 = gp(xPow2, 1.0 / (double)i); // xPow1 = xScaled^i / i! 
		
		if ((i % 2) == 1) {
			result = add(result, xPow1); 
		}
	}

    return result;
} // end of sinh()
/// <summary>sinh of bivector (uses fast special case)
/// </summary>
public static bivector sinh(bivector a)
{
	double _alpha = Math.Sqrt(Math.Abs((-a.m_e1_e2*a.m_e1_e2-a.m_e2_e3*a.m_e2_e3-a.m_e3_e1*a.m_e3_e1)));

	double _mul;
	if (_alpha != 0.0) {
		_mul = Math.Sin(_alpha)/((_alpha));

	}
	else {
		_mul = 0.0;

	}
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			_mul*a.m_e1_e2, // e1_e2
			_mul*a.m_e2_e3, // e2_e3
			_mul*a.m_e3_e1 // e3_e1
		);
}

/// <summary>Computes cosine of mv up to 12th term.
/// 
/// </summary>
public static mv cos(mv x) {
	return cos(x, 12);
}

/// <summary>Computes cosine of mv.
/// 
/// </summary>
public static mv cos(mv x, int order) {
	{ // First try special cases: check if (x * x) is scalar
		mv xSquared = gp(x, x);
		double s_xSquared = xSquared.get_scalar();
		if ((norm2_returns_scalar(xSquared) - s_xSquared * s_xSquared) < 1E-14) {
			// OK (x * x == ~scalar), so use special cases:
			if (s_xSquared > 0.0) {
				return new mv(Math.Cos(Math.Sqrt(s_xSquared)));
			}
			else if (s_xSquared < 0.0) {
				return new mv(Math.Cosh(Math.Sqrt(-s_xSquared)));
			}
			else {
				return new mv(1.0);
			}
		}
	}

	// else do general series eval . . .


	mv result = new mv(1.0);
	if (order == 0) return result;

	// taylor series approximation
	mv xPow1 = new mv(1.0); // xPow1 = 1.0
	for (int i = 1; i <= order; i++) {
		mv xPow2 = gp(xPow1, x); // xPow2 = xPow1 * x
		xPow1 = gp(xPow2, 1.0 / (double)i); // xPow1 = xScaled^i / i! 
		
		if ((i % 4) == 2)
		{
			result = subtract(result, xPow1); // result2 = result1 - xPow1
		}
		else if ((i % 4) == 0) 
		{
			result = add(result, xPow1); // result2 = result1 + xPow1
		}		
    }

	return result;
} // end of cos()
/// <summary>cos of bivector (uses fast special case)
/// </summary>
public static double cos(bivector a)
{
	double _alpha = Math.Sqrt(Math.Abs((-a.m_e1_e2*a.m_e1_e2-a.m_e2_e3*a.m_e2_e3-a.m_e3_e1*a.m_e3_e1)));

	return Math.Cosh(_alpha);
}

/// <summary>Computes sine of mv up to 12th term.
/// 
/// </summary>
public static mv sin(mv x) {
	return sin(x, 12);
}

/// <summary>Computes sine of mv.
/// 
/// </summary>
public static mv sin(mv x, int order) {
   
	{ // First try special cases: check if (x * x) is scalar
		mv xSquared = gp(x, x);
		double s_xSquared = xSquared.get_scalar();
		if ((norm2_returns_scalar(xSquared) - s_xSquared * s_xSquared) < 1E-14) {
			// OK (x * x == ~scalar), so use special cases:
			if (s_xSquared < 0.0) {
				double a = Math.Sqrt(-s_xSquared);
				return sas(x, Math.Sinh(a) / a, 0.0);
			}
			else if (s_xSquared > 0.0) {
				double a = Math.Sqrt(s_xSquared);
				return sas(x, Math.Sin(a) / a, 0.0);
			}
			else {
				return x;
			}
		}
	}

	// else do general series eval . . .

	// result = A -  ....	+ ... - ...
	mv result = new mv(); // result = 0;
    if (order == 0) return result;
    	
	// taylor series approximation
	mv xPow1 = new mv(1.0); // xPow1 = 1.0
	for (int i = 1; i <= order; i++) {
		mv xPow2 = gp(xPow1, x); // xPow2 = xPow1 * x
		
		xPow1 = gp(xPow2, 1.0 / (double)i); // xPow1 = xScaled^i / i! 
		
		if ((i % 4) == 3)
		{
			result = subtract(result, xPow1); // result = result - xPow1
		}
		else if ((i % 4) == 1) 
		{
			result = add(result, xPow1); // result = result + xPow1
		}
	}

	return result;
} // end of sin()

/// <summary>sin of bivector (uses fast special case)
/// </summary>
public static bivector sin(bivector a)
{
	double _alpha = Math.Sqrt(Math.Abs((-a.m_e1_e2*a.m_e1_e2-a.m_e2_e3*a.m_e2_e3-a.m_e3_e1*a.m_e3_e1)));

	double _mul;
	if (_alpha != 0.0) {
		_mul = Math.Sinh(_alpha)/((_alpha));

	}
	else {
		_mul = 0.0;

	}
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			_mul*a.m_e1_e2, // e1_e2
			_mul*a.m_e2_e3, // e2_e3
			_mul*a.m_e3_e1 // e3_e1
		);
}
/// <summary>Returns negation of mv.
/// </summary>
public static mv negate(mv_if a)
{
	double[][] ac = a.to_mv().c();
	double[][] cc = new double[4][];
	
	if (ac[0] != null) {
		cc[0] = new double[1];
		neg_0(ac[0], cc[0]);
	}
	
	if (ac[1] != null) {
		cc[1] = new double[3];
		neg_1(ac[1], cc[1]);
	}
	
	if (ac[2] != null) {
		cc[2] = new double[3];
		neg_2(ac[2], cc[2]);
	}
	
	if (ac[3] != null) {
		cc[3] = new double[1];
		neg_3(ac[3], cc[3]);
	}
	return new mv(cc);
}
/// <summary>Returns Clifford conjugate of mv.
/// </summary>
public static mv cliffordConjugate(mv_if a)
{
	double[][] ac = a.to_mv().c();
	double[][] cc = new double[4][];
	
	if (ac[0] != null) {
		cc[0] = new double[1];
		copyGroup_0(ac[0], cc[0]);
	}
	
	if (ac[1] != null) {
		cc[1] = new double[3];
		neg_1(ac[1], cc[1]);
	}
	
	if (ac[2] != null) {
		cc[2] = new double[3];
		neg_2(ac[2], cc[2]);
	}
	
	if (ac[3] != null) {
		cc[3] = new double[1];
		copyGroup_3(ac[3], cc[3]);
	}
	return new mv(cc);
}
/// <summary>Returns grade involution of mv.
/// </summary>
public static mv gradeInvolution(mv_if a)
{
	double[][] ac = a.to_mv().c();
	double[][] cc = new double[4][];
	
	if (ac[0] != null) {
		cc[0] = new double[1];
		copyGroup_0(ac[0], cc[0]);
	}
	
	if (ac[1] != null) {
		cc[1] = new double[3];
		neg_1(ac[1], cc[1]);
	}
	
	if (ac[2] != null) {
		cc[2] = new double[3];
		copyGroup_2(ac[2], cc[2]);
	}
	
	if (ac[3] != null) {
		cc[3] = new double[1];
		neg_3(ac[3], cc[3]);
	}
	return new mv(cc);
}
/// <summary>Returns reverse of mv.
/// </summary>
public static mv reverse(mv_if a)
{
	double[][] ac = a.to_mv().c();
	double[][] cc = new double[4][];
	
	if (ac[0] != null) {
		cc[0] = new double[1];
		copyGroup_0(ac[0], cc[0]);
	}
	
	if (ac[1] != null) {
		cc[1] = new double[3];
		copyGroup_1(ac[1], cc[1]);
	}
	
	if (ac[2] != null) {
		cc[2] = new double[3];
		neg_2(ac[2], cc[2]);
	}
	
	if (ac[3] != null) {
		cc[3] = new double[1];
		neg_3(ac[3], cc[3]);
	}
	return new mv(cc);
}
/// <summary>Returns negation of vector.
/// </summary>
public static vector negate(vector a)
{
	return new vector(vector.coord_e1_e2_e3,
			-a.m_e1, // e1
			-a.m_e2, // e2
			-a.m_e3 // e3
		);

}
/// <summary>Returns Clifford conjugate of vector.
/// </summary>
public static vector cliffordConjugate(vector a)
{
	return new vector(vector.coord_e1_e2_e3,
			-a.m_e1, // e1
			-a.m_e2, // e2
			-a.m_e3 // e3
		);

}
/// <summary>Returns grade involution of vector.
/// </summary>
public static vector gradeInvolution(vector a)
{
	return new vector(vector.coord_e1_e2_e3,
			-a.m_e1, // e1
			-a.m_e2, // e2
			-a.m_e3 // e3
		);

}
/// <summary>Returns reverse of vector.
/// </summary>
public static vector reverse(vector a)
{
	return new vector(vector.coord_e1_e2_e3,
			a.m_e1, // e1
			a.m_e2, // e2
			a.m_e3 // e3
		);

}
/// <summary>Returns negation of bivector.
/// </summary>
public static bivector negate(bivector a)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			-a.m_e1_e2, // e1_e2
			-a.m_e2_e3, // e2_e3
			-a.m_e3_e1 // e3_e1
		);

}
/// <summary>Returns Clifford conjugate of bivector.
/// </summary>
public static bivector cliffordConjugate(bivector a)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			-a.m_e1_e2, // e1_e2
			-a.m_e2_e3, // e2_e3
			-a.m_e3_e1 // e3_e1
		);

}
/// <summary>Returns grade involution of bivector.
/// </summary>
public static bivector gradeInvolution(bivector a)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			a.m_e1_e2, // e1_e2
			a.m_e2_e3, // e2_e3
			a.m_e3_e1 // e3_e1
		);

}
/// <summary>Returns reverse of bivector.
/// </summary>
public static bivector reverse(bivector a)
{
	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			-a.m_e1_e2, // e1_e2
			-a.m_e2_e3, // e2_e3
			-a.m_e3_e1 // e3_e1
		);

}
/// <summary>Returns negation of trivector.
/// </summary>
public static trivector negate(trivector a)
{
	return new trivector(trivector.coord_e1e2e3,
			-a.m_e1_e2_e3 // e1_e2_e3
		);

}
/// <summary>Returns Clifford conjugate of trivector.
/// </summary>
public static trivector cliffordConjugate(trivector a)
{
	return new trivector(trivector.coord_e1e2e3,
			a.m_e1_e2_e3 // e1_e2_e3
		);

}
/// <summary>Returns grade involution of trivector.
/// </summary>
public static trivector gradeInvolution(trivector a)
{
	return new trivector(trivector.coord_e1e2e3,
			-a.m_e1_e2_e3 // e1_e2_e3
		);

}
/// <summary>Returns reverse of trivector.
/// </summary>
public static trivector reverse(trivector a)
{
	return new trivector(trivector.coord_e1e2e3,
			-a.m_e1_e2_e3 // e1_e2_e3
		);

}
/// <summary>Returns negation of rotor.
/// </summary>
public static rotor negate(rotor a)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			-a.m_scalar, // scalar
			-a.m_e1_e2, // e1_e2
			-a.m_e2_e3, // e2_e3
			-a.m_e3_e1 // e3_e1
		);

}
/// <summary>Returns Clifford conjugate of rotor.
/// </summary>
public static rotor cliffordConjugate(rotor a)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			a.m_scalar, // scalar
			-a.m_e1_e2, // e1_e2
			-a.m_e2_e3, // e2_e3
			-a.m_e3_e1 // e3_e1
		);

}
/// <summary>Returns grade involution of rotor.
/// </summary>
public static rotor gradeInvolution(rotor a)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			a.m_scalar, // scalar
			a.m_e1_e2, // e1_e2
			a.m_e2_e3, // e2_e3
			a.m_e3_e1 // e3_e1
		);

}
/// <summary>Returns reverse of rotor.
/// </summary>
public static rotor reverse(rotor a)
{
	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			a.m_scalar, // scalar
			-a.m_e1_e2, // e1_e2
			-a.m_e2_e3, // e2_e3
			-a.m_e3_e1 // e3_e1
		);

}
/// <summary>Returns negation of e1_t.
/// </summary>
public static vector negate(e1_t a)
{
	return new vector(vector.coord_e1_e2_e3,
			-1.0, // e1
			0.0, // e2
			0.0 // e3
		);

}
/// <summary>Returns Clifford conjugate of e2_t.
/// </summary>
public static vector cliffordConjugate(e2_t a)
{
	return new vector(vector.coord_e1_e2_e3,
			0.0, // e1
			-1.0, // e2
			0.0 // e3
		);

}
/// <summary>Returns grade involution of e3_t.
/// </summary>
public static vector gradeInvolution(e3_t a)
{
	return new vector(vector.coord_e1_e2_e3,
			0.0, // e1
			0.0, // e2
			-1.0 // e3
		);

}
/// <summary>Returns reverse of I3_t.
/// </summary>
public static trivector reverse(I3_t a)
{
	return new trivector(trivector.coord_e1e2e3,
			-1.0 // e1_e2_e3
		);

}
/// <summary>Returns negation of double.
/// </summary>
public static double negate(double a)
{
	return -a;

}
/// <summary>Returns Clifford conjugate of double.
/// </summary>
public static double cliffordConjugate(double a)
{
	return a;

}
/// <summary>Returns grade involution of double.
/// </summary>
public static double gradeInvolution(double a)
{
	return a;

}
/// <summary>Returns reverse of double.
/// </summary>
public static double reverse(double a)
{
	return a;

}
/// <summary>Returns unit of mv using default metric.
/// </summary>
public static mv unit(mv_if a)
{
	double n = norm_returns_scalar(a.to_mv());
	double[][] ac = a.to_mv().c();
	double[][] cc = new double[4][];
	
	if (ac[0] != null) {
		cc[0] = new double[1];
		copyDiv_0(ac[0], cc[0], n);
	}
	
	if (ac[1] != null) {
		cc[1] = new double[3];
		copyDiv_1(ac[1], cc[1], n);
	}
	
	if (ac[2] != null) {
		cc[2] = new double[3];
		copyDiv_2(ac[2], cc[2], n);
	}
	
	if (ac[3] != null) {
		cc[3] = new double[1];
		copyDiv_3(ac[3], cc[3], n);
	}
	return new mv(cc);
}
/// <summary>Returns unit of vector using default metric.
/// </summary>
public static vector unit(vector a)
{
	double _n_ = Math.Sqrt((a.m_e1*a.m_e1+a.m_e2*a.m_e2+a.m_e3*a.m_e3));

	return new vector(vector.coord_e1_e2_e3,
			a.m_e1/((_n_)), // e1
			a.m_e2/((_n_)), // e2
			a.m_e3/((_n_)) // e3
		);
}
/// <summary>Returns unit of bivector using default metric.
/// </summary>
public static bivector unit(bivector a)
{
	double _n_ = Math.Sqrt((a.m_e1_e2*a.m_e1_e2+a.m_e2_e3*a.m_e2_e3+a.m_e3_e1*a.m_e3_e1));

	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			a.m_e1_e2/((_n_)), // e1_e2
			a.m_e2_e3/((_n_)), // e2_e3
			a.m_e3_e1/((_n_)) // e3_e1
		);
}
/// <summary>Returns unit of trivector using default metric.
/// </summary>
public static trivector unit(trivector a)
{
	double _n_ = Math.Sqrt(a.m_e1_e2_e3*a.m_e1_e2_e3);

	return new trivector(trivector.coord_e1e2e3,
			a.m_e1_e2_e3/((_n_)) // e1_e2_e3
		);
}
/// <summary>Returns unit of rotor using default metric.
/// </summary>
public static rotor unit(rotor a)
{
	double _n_ = Math.Sqrt((a.m_e1_e2*a.m_e1_e2+a.m_e2_e3*a.m_e2_e3+a.m_e3_e1*a.m_e3_e1+a.m_scalar*a.m_scalar));

	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			a.m_scalar/((_n_)), // scalar
			a.m_e1_e2/((_n_)), // e1_e2
			a.m_e2_e3/((_n_)), // e2_e3
			a.m_e3_e1/((_n_)) // e3_e1
		);
}
/// <summary>Returns unit of oddVersor using default metric.
/// </summary>
public static oddVersor unit(oddVersor a)
{
	double _n_ = Math.Sqrt((a.m_e1*a.m_e1+a.m_e1_e2_e3*a.m_e1_e2_e3+a.m_e2*a.m_e2+a.m_e3*a.m_e3));

	return new oddVersor(oddVersor.coord_e1_e2_e3_e1e2e3,
			a.m_e1/((_n_)), // e1
			a.m_e2/((_n_)), // e2
			a.m_e3/((_n_)), // e3
			a.m_e1_e2_e3/((_n_)) // e1_e2_e3
		);
}
/// <summary>Returns unit of e1_t using default metric.
/// </summary>
public static e1_t unit(e1_t a)
{
	return new e1_t(		);
}
/// <summary>Returns unit of e2_t using default metric.
/// </summary>
public static e2_t unit(e2_t a)
{
	return new e2_t(		);
}
/// <summary>Returns unit of I3_t using default metric.
/// </summary>
public static I3_t unit(I3_t a)
{
	return new I3_t(		);
}
/// <summary>Returns versor inverse of a using default metric.
/// </summary>
public static mv versorInverse(mv_if a)
{
	double n2 = norm2_returns_scalar(a.to_mv());
	double[][] ac = a.to_mv().c();
	double[][] cc = new double[4][];
	
	if (ac[0] != null) {
		cc[0] = new double[1];
		copyDiv_0(ac[0], cc[0], n2);
	}
	
	if (ac[1] != null) {
		cc[1] = new double[3];
		copyDiv_1(ac[1], cc[1], n2);
	}
	
	if (ac[2] != null) {
		cc[2] = new double[3];
		copyDiv_2(ac[2], cc[2], -n2);
	}
	
	if (ac[3] != null) {
		cc[3] = new double[1];
		copyDiv_3(ac[3], cc[3], -n2);
	}
	return new mv(cc);
}
/// <summary>Returns versor inverse of a using default metric.
/// </summary>
public static vector versorInverse(vector a)
{
	double _n2_ = (a.m_e1*a.m_e1+a.m_e2*a.m_e2+a.m_e3*a.m_e3);

	return new vector(vector.coord_e1_e2_e3,
			a.m_e1/((_n2_)), // e1
			a.m_e2/((_n2_)), // e2
			a.m_e3/((_n2_)) // e3
		);
}
/// <summary>Returns versor inverse of a using default metric.
/// </summary>
public static bivector versorInverse(bivector a)
{
	double _n2_ = (a.m_e1_e2*a.m_e1_e2+a.m_e2_e3*a.m_e2_e3+a.m_e3_e1*a.m_e3_e1);

	return new bivector(bivector.coord_e1e2_e2e3_e3e1,
			-a.m_e1_e2/((_n2_)), // e1_e2
			-a.m_e2_e3/((_n2_)), // e2_e3
			-a.m_e3_e1/((_n2_)) // e3_e1
		);
}
/// <summary>Returns versor inverse of a using default metric.
/// </summary>
public static trivector versorInverse(trivector a)
{
	double _n2_ = a.m_e1_e2_e3*a.m_e1_e2_e3;

	return new trivector(trivector.coord_e1e2e3,
			-a.m_e1_e2_e3/((_n2_)) // e1_e2_e3
		);
}
/// <summary>Returns versor inverse of a using default metric.
/// </summary>
public static rotor versorInverse(rotor a)
{
	double _n2_ = (a.m_e1_e2*a.m_e1_e2+a.m_e2_e3*a.m_e2_e3+a.m_e3_e1*a.m_e3_e1+a.m_scalar*a.m_scalar);

	return new rotor(rotor.coord_scalar_e1e2_e2e3_e3e1,
			a.m_scalar/((_n2_)), // scalar
			-a.m_e1_e2/((_n2_)), // e1_e2
			-a.m_e2_e3/((_n2_)), // e2_e3
			-a.m_e3_e1/((_n2_)) // e3_e1
		);
}
/// <summary>Returns versor inverse of a using default metric.
/// </summary>
public static e1_t versorInverse(e1_t a)
{
	return new e1_t(		);
}
/// <summary>Returns versor inverse of a using default metric.
/// </summary>
public static e3_t versorInverse(e3_t a)
{
	return new e3_t(		);
}
/// <summary>Returns versor inverse of a using default metric.
/// </summary>
public static trivector versorInverse(I3_t a)
{
	return new trivector(trivector.coord_e1e2e3,
			-1.0 // e1_e2_e3
		);
}
/// <summary>Returns true if all coordinates of a are abs <= b
/// </summary>
public static bool zero(mv_if a, double b)
{
	double[][] ac = a.to_mv().c();
	
	if (ac[0] != null) {
		if (!zeroGroup_0(ac[0], b)) return false;
	}
	
	if (ac[1] != null) {
		if (!zeroGroup_1(ac[1], b)) return false;
	}
	
	if (ac[2] != null) {
		if (!zeroGroup_2(ac[2], b)) return false;
	}
	
	if (ac[3] != null) {
		if (!zeroGroup_3(ac[3], b)) return false;
	}
	return true;
}
/// <summary>Returns true if all coordinates of a are abs <= b
/// </summary>
public static bool zero(vector a, double b)
{
	if ((a.m_e1 < -b) || (a.m_e1 > b)) return false;
	if ((a.m_e2 < -b) || (a.m_e2 > b)) return false;
	if ((a.m_e3 < -b) || (a.m_e3 > b)) return false;
	return true;
}
/// <summary>Returns true if all coordinates of a are abs <= b
/// </summary>
public static bool zero(bivector a, double b)
{
	if ((a.m_e1_e2 < -b) || (a.m_e1_e2 > b)) return false;
	if ((a.m_e2_e3 < -b) || (a.m_e2_e3 > b)) return false;
	if ((a.m_e3_e1 < -b) || (a.m_e3_e1 > b)) return false;
	return true;
}
/// <summary>Returns true if all coordinates of a are abs <= b
/// </summary>
public static bool zero(trivector a, double b)
{
	if ((a.m_e1_e2_e3 < -b) || (a.m_e1_e2_e3 > b)) return false;
	return true;
}
/// <summary>Returns true if all coordinates of a are abs <= b
/// </summary>
public static bool zero(rotor a, double b)
{
	if ((a.m_scalar < -b) || (a.m_scalar > b)) return false;
	if ((a.m_e1_e2 < -b) || (a.m_e1_e2 > b)) return false;
	if ((a.m_e2_e3 < -b) || (a.m_e2_e3 > b)) return false;
	if ((a.m_e3_e1 < -b) || (a.m_e3_e1 > b)) return false;
	return true;
}
/// <summary>Returns true if all coordinates of a are abs <= b
/// </summary>
public static bool zero(I3_t a, double b)
{
	if (1.0 > b) return false;
	return true;
}
/// <summary>Returns true if all coordinates of a are abs <= b
/// </summary>
public static bool zero(e1_t a, double b)
{
	if (1.0 > b) return false;
	return true;
}
} // end of class c3ga
} // end of namespace c3ga_ns
