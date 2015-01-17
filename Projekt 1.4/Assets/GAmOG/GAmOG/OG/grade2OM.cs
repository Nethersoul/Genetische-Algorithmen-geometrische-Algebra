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
/// Domain grade 2: e1^e2, e2^e3, -1*e1^e3.
/// Domain grade 3: .
/// 
/// The range and domain are equal.
/// 
/// </summary>
public class grade2OM 
{ 
	/// <summary>Matrix for grade 2; the size is 3 x 3
	/// </summary>
	protected internal double[] m_m2 = new double[9];

	/// <summary>Constructs a new grade2OM, set to identity.</summary>
	public grade2OM() { SetIdentity(); }

	/// <summary>Copy constructor.</summary>
	public grade2OM(grade2OM M) { Set(M); }

	/// <summary>Constructor from matrix.</summary>
	public grade2OM(double[] M) { Set(M); }

	/// <summary>Constructor from matrix.</summary>
	public grade2OM(double[] M, bool transposed) { if (transposed) { SetTranspose(M); } else { Set(M); } }
	
	/// <summary>Constructor from images of basis vectors.</summary>
	public grade2OM(vector ie1, vector ie2, vector ie3)
		{ Set(ie1, ie2, ie3); }

	/// <summary>Converts a om to a grade2OM.
	/// Warning: coordinates which cannot be represented are silenty lost.</summary>
	public grade2OM(om M) { Set(M); }

/// <summary>Sets grade2OM to identity.
/// </summary>
public void SetIdentity()
{
	m_m2[0] = 1.0;
	m_m2[3] = m_m2[6] = 0.0;

	m_m2[1] = m_m2[7] = 0.0;
	m_m2[4] = 1.0;

	m_m2[2] = m_m2[5] = 0.0;
	m_m2[8] = 1.0;

}

/// <summary>Copies grade2OM.
/// </summary>
public void Set(grade2OM src)
{
	m_m2[0] = src.m_m2[0];
	m_m2[3] = src.m_m2[3];
	m_m2[6] = src.m_m2[6];

	m_m2[1] = src.m_m2[1];
	m_m2[4] = src.m_m2[4];
	m_m2[7] = src.m_m2[7];

	m_m2[2] = src.m_m2[2];
	m_m2[5] = src.m_m2[5];
	m_m2[8] = src.m_m2[8];

}
/// <summary>Sets grade2OM from images of the domain vectors.
/// </summary>
public void Set(vector ie1, vector ie2, vector ie3)
{
	m_m2[0] = (ie1.m_e1*ie2.m_e2-ie1.m_e2*ie2.m_e1);
	m_m2[3] = (ie1.m_e2*ie2.m_e3-ie1.m_e3*ie2.m_e2);
	m_m2[6] = -(ie1.m_e1*ie2.m_e3-ie1.m_e3*ie2.m_e1);

	m_m2[1] = (ie2.m_e1*ie3.m_e2-ie2.m_e2*ie3.m_e1);
	m_m2[4] = (ie2.m_e2*ie3.m_e3-ie2.m_e3*ie3.m_e2);
	m_m2[7] = -(ie2.m_e1*ie3.m_e3-ie2.m_e3*ie3.m_e1);

	m_m2[2] = (-ie1.m_e1*ie3.m_e2+ie1.m_e2*ie3.m_e1);
	m_m2[5] = (-ie1.m_e2*ie3.m_e3+ie1.m_e3*ie3.m_e2);
	m_m2[8] = -(-ie1.m_e1*ie3.m_e3+ie1.m_e3*ie3.m_e1);

}
/// <summary>Sets grade2OM from a matrix.
/// </summary>
public void Set(double[] M)
{
	m_m2[0] = (M[0]*M[4]-M[1]*M[3]);
	m_m2[3] = (M[3]*M[7]-M[4]*M[6]);
	m_m2[6] = -(M[0]*M[7]-M[1]*M[6]);

	m_m2[1] = (M[1]*M[5]-M[2]*M[4]);
	m_m2[4] = (M[4]*M[8]-M[5]*M[7]);
	m_m2[7] = -(M[1]*M[8]-M[2]*M[7]);

	m_m2[2] = (-M[0]*M[5]+M[2]*M[3]);
	m_m2[5] = (-M[3]*M[8]+M[5]*M[6]);
	m_m2[8] = -(-M[0]*M[8]+M[2]*M[6]);

}
/// <summary>Sets grade2OM from a transposed matrix.
/// </summary>
public void SetTranspose(double[] M)
{
	m_m2[0] = (M[0]*M[4]-M[1]*M[3]);
	m_m2[3] = (M[1]*M[5]-M[2]*M[4]);
	m_m2[6] = -(M[0]*M[5]-M[2]*M[3]);

	m_m2[1] = (M[3]*M[7]-M[4]*M[6]);
	m_m2[4] = (M[4]*M[8]-M[5]*M[7]);
	m_m2[7] = -(M[3]*M[8]-M[5]*M[6]);

	m_m2[2] = (-M[0]*M[7]+M[1]*M[6]);
	m_m2[5] = (-M[1]*M[8]+M[2]*M[7]);
	m_m2[8] = -(-M[0]*M[8]+M[2]*M[6]);

}
/// <summary>Copies a om to a grade2OM
/// Warning 1: coordinates which cannot be represented are silenty lost.
/// Warning 2: coordinates which are not present in 'src' are set to zero in 'dst'.
/// 
/// </summary>
void Set(om src) {
	m_m2[0] =  src.m_m2[0];
	m_m2[1] =  src.m_m2[2];
	m_m2[2] = -1.0 *  src.m_m2[1];
	m_m2[3] =  src.m_m2[6];
	m_m2[4] =  src.m_m2[8];
	m_m2[5] = -1.0 *  src.m_m2[7];
	m_m2[6] = -1.0 *  src.m_m2[3];
	m_m2[7] = -1.0 *  src.m_m2[5];
	m_m2[8] =  src.m_m2[4];
}
} // end of class grade2OM
} // end of namespace c3ga_ns
