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
/// <summary>This class can hold a specialized outermorphism.
/// 
/// The coordinates are stored in type double.
/// 
/// There are 4 matrices, one for each grade.
/// The columns of these matrices are the range of the outermorphism.
/// Matrices are stored in row-major order. So the coordinates of rows are stored contiguously.
/// Domain grade 1: e1, e2, e3.
/// Domain grade 2: .
/// Domain grade 3: .
/// 
/// The range and domain are equal.
/// 
/// </summary>
public class grade1OM 
{ 
	/// <summary>Matrix for grade 1; the size is 3 x 3
	/// </summary>
	protected internal double[] m_m1 = new double[9];

	/// <summary>Constructs a new grade1OM, set to identity.</summary>
	public grade1OM() { SetIdentity(); }

	/// <summary>Copy constructor.</summary>
	public grade1OM(grade1OM M) { Set(M); }

	/// <summary>Constructor from matrix.</summary>
	public grade1OM(double[] M) { Set(M); }

	/// <summary>Constructor from matrix.</summary>
	public grade1OM(double[] M, bool transposed) { if (transposed) { SetTranspose(M); } else { Set(M); } }
	
	/// <summary>Constructor from images of basis vectors.</summary>
	public grade1OM(vector ie1, vector ie2, vector ie3)
		{ Set(ie1, ie2, ie3); }

	/// <summary>Converts a om to a grade1OM.
	/// Warning: coordinates which cannot be represented are silenty lost.</summary>
	public grade1OM(om M) { Set(M); }

/// <summary>Sets grade1OM to identity.
/// </summary>
public void SetIdentity()
{
	m_m1[0] = 1.0;
	m_m1[3] = m_m1[6] = 0.0;

	m_m1[1] = m_m1[7] = 0.0;
	m_m1[4] = 1.0;

	m_m1[2] = m_m1[5] = 0.0;
	m_m1[8] = 1.0;

}

/// <summary>Copies grade1OM.
/// </summary>
public void Set(grade1OM src)
{
	m_m1[0] = src.m_m1[0];
	m_m1[3] = src.m_m1[3];
	m_m1[6] = src.m_m1[6];

	m_m1[1] = src.m_m1[1];
	m_m1[4] = src.m_m1[4];
	m_m1[7] = src.m_m1[7];

	m_m1[2] = src.m_m1[2];
	m_m1[5] = src.m_m1[5];
	m_m1[8] = src.m_m1[8];

}
/// <summary>Sets grade1OM from images of the domain vectors.
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

}
/// <summary>Sets grade1OM from a matrix.
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

}
/// <summary>Sets grade1OM from a transposed matrix.
/// </summary>
public void SetTranspose(double[] M)
{
	m_m1[0] = M[0];
	m_m1[3] = M[1];
	m_m1[6] = M[2];

	m_m1[1] = M[3];
	m_m1[4] = M[4];
	m_m1[7] = M[5];

	m_m1[2] = M[6];
	m_m1[5] = M[7];
	m_m1[8] = M[8];

}
/// <summary>Copies a om to a grade1OM
/// Warning 1: coordinates which cannot be represented are silenty lost.
/// Warning 2: coordinates which are not present in 'src' are set to zero in 'dst'.
/// 
/// </summary>
void Set(om src) {
	m_m1[0] =  src.m_m1[0];
	m_m1[1] =  src.m_m1[1];
	m_m1[2] =  src.m_m1[2];
	m_m1[3] =  src.m_m1[3];
	m_m1[4] =  src.m_m1[4];
	m_m1[5] =  src.m_m1[5];
	m_m1[6] =  src.m_m1[6];
	m_m1[7] =  src.m_m1[7];
	m_m1[8] =  src.m_m1[8];
}
} // end of class grade1OM
} // end of namespace c3ga_ns
