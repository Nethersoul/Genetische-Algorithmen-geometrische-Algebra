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
/// <summary>This class can hold a general multivector.
/// 
/// The coordinates are stored in type double.
/// 
/// There are 4 coordinate groups:
/// group 0:1  (grade 0).
/// group 1:e1, e2, e3  (grade 1).
/// group 2:e1^e2, e1^e3, e2^e3  (grade 2).
/// group 3:e1^e2^e3  (grade 3).
/// 
/// 8 doubles are allocated inside the struct.
/// 
/// </summary>
public class mv  :  mv_if
{ 

    /// <summary>
	/// the coordinates
    /// </summary>
	protected internal double[][] m_c = new double[4][]; 
	


    /// <summary>
	/// Constructs a new mv with value 0.
    /// </summary>
	public mv() {Set();}

    /// <summary>
	/// Copy constructor.
    /// </summary>
	public mv(mv A) {Set(A);}


    /// <summary>
	/// Constructs a new mv with scalar value 'scalar'.
    /// </summary>
	public mv(double scalar) {Set(scalar);}

    /// <summary>
	/// Constructs a new mv from compressed 'coordinates'.
	/// <param name="gu">bitwise OR of the GRADEs or GROUPs that are non-zero.</param>
	/// <param name="coordinates"> compressed coordinates.</param>
    /// </summary>
	public mv(GroupBitmap gu, double[] coordinates) {Set(gu, coordinates);}

    /// <summary>
	/// Constructs a new mv from 'coordinates'.
	/// <param name="coordinates">The coordinates (one array for each group, entries may be null). The arrays are kept.</param>
    /// </summary>
	public mv(double[][] coordinates) {Set(coordinates);}
	
    /// <summary>
	/// Converts a e1_t to a mv.
    /// </summary>
	public mv(e1_t A) {Set(A);}
    /// <summary>
	/// Converts a e2_t to a mv.
    /// </summary>
	public mv(e2_t A) {Set(A);}
    /// <summary>
	/// Converts a e3_t to a mv.
    /// </summary>
	public mv(e3_t A) {Set(A);}
    /// <summary>
	/// Converts a I3_t to a mv.
    /// </summary>
	public mv(I3_t A) {Set(A);}
    /// <summary>
	/// Converts a vector to a mv.
    /// </summary>
	public mv(vector A) {Set(A);}
    /// <summary>
	/// Converts a bivector to a mv.
    /// </summary>
	public mv(bivector A) {Set(A);}
    /// <summary>
	/// Converts a trivector to a mv.
    /// </summary>
	public mv(trivector A) {Set(A);}
    /// <summary>
	/// Converts a rotor to a mv.
    /// </summary>
	public mv(rotor A) {Set(A);}
    /// <summary>
	/// Converts a oddVersor to a mv.
    /// </summary>
	public mv(oddVersor A) {Set(A);}


    /// <summary>
	/// returns group usage bitmap
    /// </summary>
	public GroupBitmap gu() {
		return 
			((m_c[0] == null) ? 0 : GroupBitmap.GROUP_0) |
			((m_c[1] == null) ? 0 : GroupBitmap.GROUP_1) |
			((m_c[2] == null) ? 0 : GroupBitmap.GROUP_2) |
			((m_c[3] == null) ? 0 : GroupBitmap.GROUP_3) |
			0;
	}
    /// <summary>
	/// Returns array of array of coordinates.
	/// Each entry contain the coordinates for one group/grade.
    /// </summary>
	public double[][] c() { return m_c; }
	
	/// <summary>sets this to 0.
	/// </summary>
	public void Set() {
		m_c[0] = null;
		m_c[1] = null;
		m_c[2] = null;
		m_c[3] = null;
	}
	/// <summary>sets this to scalar value.
	/// </summary>
	public void Set(double val) {
		AllocateGroups(GroupBitmap.GROUP_0);
		m_c[0][0] = val;
	}
	/// <summary>sets this coordinates in 'arr'.
	/// </summary>
	/// <param name="gu">bitwise or of the GROUPs and GRADEs which are present in 'arr'.
	/// </param>
	/// <param name="arr">compressed coordinates.
	/// </param>
	public void Set(GroupBitmap gu, double[] arr) {
		AllocateGroups(gu);
		int idx = 0;
		if ((gu & GroupBitmap.GROUP_0) != 0) {
			for (int i = 0; i < 1; i++)
				m_c[0][i] = arr[idx + i];
			idx += 1;
		}
		if ((gu & GroupBitmap.GROUP_1) != 0) {
			for (int i = 0; i < 3; i++)
				m_c[1][i] = arr[idx + i];
			idx += 3;
		}
		if ((gu & GroupBitmap.GROUP_2) != 0) {
			for (int i = 0; i < 3; i++)
				m_c[2][i] = arr[idx + i];
			idx += 3;
		}
		if ((gu & GroupBitmap.GROUP_3) != 0) {
			for (int i = 0; i < 1; i++)
				m_c[3][i] = arr[idx + i];
			idx += 1;
		}
	}
	/// <summary>sets this coordinates in 'arr'. 
	/// 'arr' is kept, so changes to 'arr' will be reflected in the value of this multivector. Make sure 'arr' has length 4 and each subarray has the length of the respective group/grade
	/// </summary>
	/// <param name="arr">coordinates.
	/// </param>
	public void Set(double[][] arr) {
		m_c = arr;
	}
	/// <summary>sets this to multivector value.
	/// </summary>
	public void Set(mv src) {
		AllocateGroups(src.gu());
		if (m_c[0] != null) {
			c3ga.Copy_1(m_c[0], src.m_c[0]);
		}
		if (m_c[1] != null) {
			c3ga.Copy_3(m_c[1], src.m_c[1]);
		}
		if (m_c[2] != null) {
			c3ga.Copy_3(m_c[2], src.m_c[2]);
		}
		if (m_c[3] != null) {
			c3ga.Copy_1(m_c[3], src.m_c[3]);
		}
	}

	/// <summary>sets this to e1_t value.
	/// </summary>
	public void Set(e1_t src) {
		AllocateGroups(GroupBitmap.GROUP_1);
		double[] ptr;

		ptr = m_c[1];
		ptr[0] = 1.0;
		ptr[1] = ptr[2] = 0.0;
	}

	/// <summary>sets this to e2_t value.
	/// </summary>
	public void Set(e2_t src) {
		AllocateGroups(GroupBitmap.GROUP_1);
		double[] ptr;

		ptr = m_c[1];
		ptr[0] = ptr[2] = 0.0;
		ptr[1] = 1.0;
	}

	/// <summary>sets this to e3_t value.
	/// </summary>
	public void Set(e3_t src) {
		AllocateGroups(GroupBitmap.GROUP_1);
		double[] ptr;

		ptr = m_c[1];
		ptr[0] = ptr[1] = 0.0;
		ptr[2] = 1.0;
	}

	/// <summary>sets this to I3_t value.
	/// </summary>
	public void Set(I3_t src) {
		AllocateGroups(GroupBitmap.GROUP_3);
		double[] ptr;

		ptr = m_c[3];
		ptr[0] = 1.0;
	}

	/// <summary>sets this to vector value.
	/// </summary>
	public void Set(vector src) {
		AllocateGroups(GroupBitmap.GROUP_1);
		double[] ptr;

		ptr = m_c[1];
		ptr[0] = src.m_e1;
		ptr[1] = src.m_e2;
		ptr[2] = src.m_e3;
	}

	/// <summary>sets this to bivector value.
	/// </summary>
	public void Set(bivector src) {
		AllocateGroups(GroupBitmap.GROUP_2);
		double[] ptr;

		ptr = m_c[2];
		ptr[0] = src.m_e1_e2;
		ptr[1] = -src.m_e3_e1;
		ptr[2] = src.m_e2_e3;
	}

	/// <summary>sets this to trivector value.
	/// </summary>
	public void Set(trivector src) {
		AllocateGroups(GroupBitmap.GROUP_3);
		double[] ptr;

		ptr = m_c[3];
		ptr[0] = src.m_e1_e2_e3;
	}

	/// <summary>sets this to rotor value.
	/// </summary>
	public void Set(rotor src) {
		AllocateGroups(GroupBitmap.GROUP_0|GroupBitmap.GROUP_2);
		double[] ptr;

		ptr = m_c[0];
		ptr[0] = src.m_scalar;

		ptr = m_c[2];
		ptr[0] = src.m_e1_e2;
		ptr[1] = -src.m_e3_e1;
		ptr[2] = src.m_e2_e3;
	}

	/// <summary>sets this to oddVersor value.
	/// </summary>
	public void Set(oddVersor src) {
		AllocateGroups(GroupBitmap.GROUP_1|GroupBitmap.GROUP_3);
		double[] ptr;

		ptr = m_c[1];
		ptr[0] = src.m_e1;
		ptr[1] = src.m_e2;
		ptr[2] = src.m_e3;

		ptr = m_c[3];
		ptr[0] = src.m_e1_e2_e3;
	}
	/// <summary>Returns the scalar coordinate of this mv
	/// </summary>
	public double get_scalar()  {
		return (m_c[0] == null) ? 0.0: m_c[0][0];
	}
	/// <summary>Returns the e1 coordinate of this mv
	/// </summary>
	public double get_e1()  {
		return (m_c[1] == null) ? 0.0: m_c[1][0];
	}
	/// <summary>Returns the e2 coordinate of this mv
	/// </summary>
	public double get_e2()  {
		return (m_c[1] == null) ? 0.0: m_c[1][1];
	}
	/// <summary>Returns the e3 coordinate of this mv
	/// </summary>
	public double get_e3()  {
		return (m_c[1] == null) ? 0.0: m_c[1][2];
	}
	/// <summary>Returns the e1_e2 coordinate of this mv
	/// </summary>
	public double get_e1_e2()  {
		return (m_c[2] == null) ? 0.0: m_c[2][0];
	}
	/// <summary>Returns the e1_e3 coordinate of this mv
	/// </summary>
	public double get_e1_e3()  {
		return (m_c[2] == null) ? 0.0: m_c[2][1];
	}
	/// <summary>Returns the e2_e3 coordinate of this mv
	/// </summary>
	public double get_e2_e3()  {
		return (m_c[2] == null) ? 0.0: m_c[2][2];
	}
	/// <summary>Returns the e1_e2_e3 coordinate of this mv
	/// </summary>
	public double get_e1_e2_e3()  {
		return (m_c[3] == null) ? 0.0: m_c[3][0];
	}

    /// <summary>
	/// Reserves memory for the groups specified by 'gu'.
	/// Keeps old memory (and values) when possible.
    /// </summary>
	private void AllocateGroups(GroupBitmap gu) {
		for (int i = 0; (1 << i) <= (int)gu; i++) {
			if (((1 << i) & (int)gu) != 0) {
				if (m_c[i] == null)
					m_c[i] = new double[c3ga.MvSize[1 << i]];
			}
			else m_c[i] = null;
		}		
	}

	/// <summary>
	/// Reserves memory for coordinate GROUP_0.
	/// If the group is already present, nothing changes.
	/// If the group is not present, memory is allocated for the new group,
	/// and the coordinates for the group are set to zero.
	/// </summary>
	private void ReserveGroup_0() {
		if (m_c[0] == null) {
			m_c[0] = new double[1];
		}
	}
	/// <summary>
	/// Reserves memory for coordinate GROUP_1.
	/// If the group is already present, nothing changes.
	/// If the group is not present, memory is allocated for the new group,
	/// and the coordinates for the group are set to zero.
	/// </summary>
	private void ReserveGroup_1() {
		if (m_c[1] == null) {
			m_c[1] = new double[3];
		}
	}
	/// <summary>
	/// Reserves memory for coordinate GROUP_2.
	/// If the group is already present, nothing changes.
	/// If the group is not present, memory is allocated for the new group,
	/// and the coordinates for the group are set to zero.
	/// </summary>
	private void ReserveGroup_2() {
		if (m_c[2] == null) {
			m_c[2] = new double[3];
		}
	}
	/// <summary>
	/// Reserves memory for coordinate GROUP_3.
	/// If the group is already present, nothing changes.
	/// If the group is not present, memory is allocated for the new group,
	/// and the coordinates for the group are set to zero.
	/// </summary>
	private void ReserveGroup_3() {
		if (m_c[3] == null) {
			m_c[3] = new double[1];
		}
	}
	/// Sets the scalar coordinate of this mv.
	public void set_scalar(double val)  {
		ReserveGroup_0();
		m_c[0][0] =  val;
	}
	/// Sets the e1 coordinate of this mv.
	public void set_e1(double val)  {
		ReserveGroup_1();
		m_c[1][0] =  val;
	}
	/// Sets the e2 coordinate of this mv.
	public void set_e2(double val)  {
		ReserveGroup_1();
		m_c[1][1] =  val;
	}
	/// Sets the e3 coordinate of this mv.
	public void set_e3(double val)  {
		ReserveGroup_1();
		m_c[1][2] =  val;
	}
	/// Sets the e1_e2 coordinate of this mv.
	public void set_e1_e2(double val)  {
		ReserveGroup_2();
		m_c[2][0] =  val;
	}
	/// Sets the e1_e3 coordinate of this mv.
	public void set_e1_e3(double val)  {
		ReserveGroup_2();
		m_c[2][1] =  val;
	}
	/// Sets the e2_e3 coordinate of this mv.
	public void set_e2_e3(double val)  {
		ReserveGroup_2();
		m_c[2][2] =  val;
	}
	/// Sets the e1_e2_e3 coordinate of this mv.
	public void set_e1_e2_e3(double val)  {
		ReserveGroup_3();
		m_c[3][0] =  val;
	}

	/// <summary>returns the absolute largest coordinate.</summary>
	public double LargestCoordinate() {
		double maxValue = 0.0, C;
		for (int g = 0; g < m_c.Length; g++) {
			if (m_c[g] != null) {
				double[] Cg = m_c[g];
				for (int b = 0; b < Cg.Length; b++) {
					C = Math.Abs(Cg[b]);
					if (C > maxValue) {
						maxValue = C;
					}
				}
			}
		}
		return maxValue;
	}
	
	/// <summary>returns the absolute largest coordinate and the corresponding basis blade bitmap (in 'bm') .</summary>
	public double LargestBasisBlade(ref int bm) {
		double maxC = -1.0, C;

		int idx = 0; // global index into coordinates (run from 0 to 8).
		bm = 0;
		
		for (int g = 0; g < m_c.Length; g++) {
			if (m_c[g] != null) {
				double[] Cg = m_c[g];
				for (int b = 0; b < m_c[g].Length; b++) {
					C = Math.Abs(Cg[b]);
					if (C > maxC) {
						maxC = C;
						bm = c3ga.BasisElementBitmapByIndex[idx];
					}
					idx++;
				}
			
			}
			else idx += c3ga.GroupSize[g];
		}

		return maxC;
	} // end of LargestBasisBlade()
	
	/// <summary>Releases memory for (near-)zero groups/grades.
	/// This also speeds up subsequent operations, because those do not have to process the released groups/grades anymore.
	/// </summary>
	/// <param name="eps">A positive threshold value.
	/// Coordinates which are smaller than epsilon are considered to be zero.
	/// </param>
	public void Compress(double eps)  {
		if ((m_c[0] != null) && c3ga.zeroGroup_0(m_c[0], eps))
			m_c[0] = null;
		if ((m_c[1] != null) && c3ga.zeroGroup_1(m_c[1], eps))
			m_c[1] = null;
		if ((m_c[2] != null) && c3ga.zeroGroup_2(m_c[2], eps))
			m_c[2] = null;
		if ((m_c[3] != null) && c3ga.zeroGroup_3(m_c[3], eps))
			m_c[3] = null;
	}

	/// <summary>
	/// Returns this multivector, converted to a string.
	/// The floating point formatter is controlled via c3ga.setStringFormat().
	/// </summary>
	public override string ToString() {
		return c3ga.String(this);
	}
	
	/// <summary>
	/// Returns this multivector, converted to a string.
	/// The floating point formatter is "F".
	/// </summary>
	public string ToString_f() {
		return ToString("F");
	}
	
	/// <summary>
	/// Returns this multivector, converted to a string.
	/// The floating point formatter is "E".
	/// </summary>
	public string ToString_e() {
		return ToString("E");
	}
	
	/// <summary>
	/// Returns this multivector, converted to a string.
	/// The floating point formatter is "E20".
	/// </summary>
	public string ToString_e20() {
		return ToString("E20");
	}
	
	/// <summary>
	/// Returns this multivector, converted to a string.
	/// <param name="fp">floating point format. Use 'null' for the default format (see setStringFormat()).</param>
	/// </summary>
	public string ToString(string fp) {
		return c3ga.String(this, fp);
	}

    /// <summary>	
    /// Converts this multivector to a 'mv' (implementation of interface 'mv_interface')
    /// </summary>
    public mv to_mv()
    {
        return this;
    }
} // end of class mv
} // end of namespace c3ga_ns
