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
namespace c3ga_ns {
/// <summary>This class can hold a general outermorphism.
/// 
/// The coordinates are stored in type double.
/// 
/// There are 4 matrices, one for each grade.
/// The columns of these matrices are the range of the outermorphism.
/// Matrices are stored in row-major order. So the coordinates of rows are stored contiguously.
/// Domain grade 1: e1, e2, e3.
/// Domain grade 2: e1^e2, e1^e3, e2^e3.
/// Domain grade 3: e1^e2^e3.
/// 
/// The range and domain are equal.
/// 
/// </summary>
public class om 
{ 
	/// <summary>Matrix for grade 1; the size is 3 x 3
	/// </summary>
	protected internal double[] m_m1 = new double[9];
	/// <summary>Matrix for grade 2; the size is 3 x 3
	/// </summary>
	protected internal double[] m_m2 = new double[9];
	/// <summary>Matrix for grade 3; the size is 1 x 1
	/// </summary>
	protected internal double[] m_m3 = new double[1];

	/// <summary>Constructs a new om, set to identity.</summary>
	public om() { SetIdentity(); }

	/// <summary>Copy constructor.</summary>
	public om(om M) { Set(M); }

	/// <summary>Constructor from matrix.</summary>
	public om(double[] M) { Set(M); }

	/// <summary>Constructor from images of basis vectors.</summary>
	public om( vector ie1,  vector ie2,  vector ie3)
		{ Set(ie1, ie2, ie3); }

	/// <summary>Converts a grade1OM to a om.
	/// Warning 1: coordinates which cannot be represented are silenty lost.
	/// Warning 2: coordinates which are not present in 'src' are set to zero in 'dst'.</summary>
	public om(grade1OM M) { Set(M); }
	
	/// <summary>Converts a grade2OM to a om.
	/// Warning 1: coordinates which cannot be represented are silenty lost.
	/// Warning 2: coordinates which are not present in 'src' are set to zero in 'dst'.</summary>
	public om(grade2OM M) { Set(M); }
	
	/// <summary>Converts a grade3OM to a om.
	/// Warning 1: coordinates which cannot be represented are silenty lost.
	/// Warning 2: coordinates which are not present in 'src' are set to zero in 'dst'.</summary>
	public om(grade3OM M) { Set(M); }
	

	public void SetIdentity() {
		c3ga.Zero_9(m_m1);

		c3ga.Zero_9(m_m2);

		c3ga.Zero_1(m_m3);

		m_m1[0] = m_m1[4] = m_m1[8] = m_m2[0] = m_m2[4] = m_m2[8] = m_m3[0] = 1.0;
	}

	void Set(om src) {
		c3ga.Copy_9(m_m1, src.m_m1);

		c3ga.Copy_9(m_m2, src.m_m2);

		c3ga.Copy_1(m_m3, src.m_m3);

	}
/// <summary>Sets grade 2 part of outermorphism matrix based on lower grade parts.
/// </summary>
public void Set_grade_2_0()
{
	m_m2[0] = (m_m1[0]*m_m1[4]-m_m1[1]*m_m1[3]);
	m_m2[3] = (m_m1[0]*m_m1[7]-m_m1[1]*m_m1[6]);
	m_m2[6] = (m_m1[3]*m_m1[7]-m_m1[4]*m_m1[6]);

}
/// <summary>Sets grade 2 part of outermorphism matrix based on lower grade parts.
/// </summary>
public void Set_grade_2_1()
{
	m_m2[1] = (m_m1[0]*m_m1[5]-m_m1[2]*m_m1[3]);
	m_m2[4] = (m_m1[0]*m_m1[8]-m_m1[2]*m_m1[6]);
	m_m2[7] = (m_m1[3]*m_m1[8]-m_m1[5]*m_m1[6]);

}
/// <summary>Sets grade 2 part of outermorphism matrix based on lower grade parts.
/// </summary>
public void Set_grade_2_2()
{
	m_m2[2] = (m_m1[1]*m_m1[5]-m_m1[2]*m_m1[4]);
	m_m2[5] = (m_m1[1]*m_m1[8]-m_m1[2]*m_m1[7]);
	m_m2[8] = (m_m1[4]*m_m1[8]-m_m1[5]*m_m1[7]);

}
/// <summary>Sets grade 3 part of outermorphism matrix based on lower grade parts.
/// </summary>
public void Set_grade_3_0()
{
	m_m3[0] = (m_m1[0]*m_m2[8]-m_m1[3]*m_m2[5]+m_m1[6]*m_m2[2]);

}
/// <summary>Sets om from images of the domain vectors.
/// </summary>
public void Set(vector ie1, vector ie2, vector ie3)
{
	m_m1[0] = ie1.m_e1;
	m_m1[3] = ie1.m_e2;
	m_m1[6] = ie1.m_e3;

	m_m1[1] = ie2.m_e1;
	m_m1[4] = ie2.m_e2;
	m_m1[7] = ie2.m_e3;

	m_m1[2] = ie3.m_e1;
	m_m1[5] = ie3.m_e2;
	m_m1[8] = ie3.m_e3;

	Set_grade_2_0();
	Set_grade_2_1();
	Set_grade_2_2();
	Set_grade_3_0();
}
/// <summary>Sets om from a matrix
/// </summary>
public void Set(double[] M)
{
	m_m1[0] = M[0];
	m_m1[3] = M[3];
	m_m1[6] = M[6];

	m_m1[1] = M[1];
	m_m1[4] = M[4];
	m_m1[7] = M[7];

	m_m1[2] = M[2];
	m_m1[5] = M[5];
	m_m1[8] = M[8];

	Set_grade_2_0();
	Set_grade_2_1();
	Set_grade_2_2();
	Set_grade_3_0();
}
/// <summary>Copies a grade1OM to a om
/// Warning 1: coordinates which cannot be represented are silenty lost.
/// Warning 2: coordinates which are not present in 'src' are set to zero in 'dst'.
/// 
/// </summary>
void Set(grade1OM src) {
	m_m1[0] =  src.m_m1[0];
	m_m1[1] =  src.m_m1[1];
	m_m1[2] =  src.m_m1[2];
	m_m1[3] =  src.m_m1[3];
	m_m1[4] =  src.m_m1[4];
	m_m1[5] =  src.m_m1[5];
	m_m1[6] =  src.m_m1[6];
	m_m1[7] =  src.m_m1[7];
	m_m1[8] =  src.m_m1[8];
	m_m2[0] = m_m2[1] = m_m2[2] = m_m2[3] = m_m2[4] = m_m2[5] = m_m2[6] = m_m2[7] = m_m2[8] = 
		m_m3[0] = 0.0;
}
/// <summary>Copies a grade2OM to a om
/// Warning 1: coordinates which cannot be represented are silenty lost.
/// Warning 2: coordinates which are not present in 'src' are set to zero in 'dst'.
/// 
/// </summary>
void Set(grade2OM src) {
	m_m2[0] =  src.m_m2[0];
	m_m2[1] = -1.0 *  src.m_m2[2];
	m_m2[2] =  src.m_m2[1];
	m_m2[3] = -1.0 *  src.m_m2[6];
	m_m2[4] =  src.m_m2[8];
	m_m2[5] = -1.0 *  src.m_m2[7];
	m_m2[6] =  src.m_m2[3];
	m_m2[7] = -1.0 *  src.m_m2[5];
	m_m2[8] =  src.m_m2[4];
	m_m1[0] = m_m1[1] = m_m1[2] = m_m1[3] = m_m1[4] = m_m1[5] = m_m1[6] = m_m1[7] = m_m1[8] = 
		m_m3[0] = 0.0;
}
/// <summary>Copies a grade3OM to a om
/// Warning 1: coordinates which cannot be represented are silenty lost.
/// Warning 2: coordinates which are not present in 'src' are set to zero in 'dst'.
/// 
/// </summary>
void Set(grade3OM src) {
	m_m3[0] =  src.m_m3[0];
	m_m1[0] = m_m1[1] = m_m1[2] = m_m1[3] = m_m1[4] = m_m1[5] = m_m1[6] = m_m1[7] = m_m1[8] = 
		m_m2[0] = m_m2[1] = m_m2[2] = m_m2[3] = m_m2[4] = m_m2[5] = m_m2[6] = m_m2[7] = m_m2[8] = 
		0.0;
}
} // end of class om
} // end of namespace c3ga_ns
