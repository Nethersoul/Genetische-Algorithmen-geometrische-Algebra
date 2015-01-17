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
/// Domain grade 1: .
/// Domain grade 2: .
/// Domain grade 3: e1^e2^e3.
/// 
/// The range and domain are equal.
/// 
/// </summary>
public class grade3OM 
{ 
	/// <summary>Matrix for grade 3; the size is 1 x 1
	/// </summary>
	protected internal double[] m_m3 = new double[1];

	/// <summary>Constructs a new grade3OM, set to identity.</summary>
	public grade3OM() { SetIdentity(); }

	/// <summary>Copy constructor.</summary>
	public grade3OM(grade3OM M) { Set(M); }

	/// <summary>Constructor from matrix.</summary>
	public grade3OM(double[] M) { Set(M); }

	/// <summary>Constructor from matrix.</summary>
	public grade3OM(double[] M, bool transposed) { if (transposed) { SetTranspose(M); } else { Set(M); } }
	
	/// <summary>Constructor from images of basis vectors.</summary>
	public grade3OM(vector ie1, vector ie2, vector ie3)
		{ Set(ie1, ie2, ie3); }

	/// <summary>Converts a om to a grade3OM.
	/// Warning: coordinates which cannot be represented are silenty lost.</summary>
	public grade3OM(om M) { Set(M); }

/// <summary>Sets grade3OM to identity.
/// </summary>
public void SetIdentity()
{
	m_m3[0] = 1.0;

}

/// <summary>Copies grade3OM.
/// </summary>
public void Set(grade3OM src)
{
	m_m3[0] = src.m_m3[0];

}
/// <summary>Sets grade3OM from images of the domain vectors.
/// </summary>
public void Set(vector ie1, vector ie2, vector ie3)
{
	m_m3[0] = (ie1.m_e1*ie2.m_e2*ie3.m_e3-ie1.m_e1*ie2.m_e3*ie3.m_e2-ie1.m_e2*ie2.m_e1*ie3.m_e3+ie1.m_e2*ie2.m_e3*ie3.m_e1+ie1.m_e3*ie2.m_e1*ie3.m_e2-ie1.m_e3*ie2.m_e2*ie3.m_e1);

}
/// <summary>Sets grade3OM from a matrix.
/// </summary>
public void Set(double[] M)
{
	m_m3[0] = (M[0]*M[4]*M[8]-M[0]*M[5]*M[7]-M[1]*M[3]*M[8]+M[1]*M[5]*M[6]+M[2]*M[3]*M[7]-M[2]*M[4]*M[6]);

}
/// <summary>Sets grade3OM from a transposed matrix.
/// </summary>
public void SetTranspose(double[] M)
{
	m_m3[0] = (M[0]*M[4]*M[8]-M[0]*M[5]*M[7]-M[1]*M[3]*M[8]+M[1]*M[5]*M[6]+M[2]*M[3]*M[7]-M[2]*M[4]*M[6]);

}
/// <summary>Copies a om to a grade3OM
/// Warning 1: coordinates which cannot be represented are silenty lost.
/// Warning 2: coordinates which are not present in 'src' are set to zero in 'dst'.
/// 
/// </summary>
void Set(om src) {
	m_m3[0] =  src.m_m3[0];
}
} // end of class grade3OM
} // end of namespace c3ga_ns
