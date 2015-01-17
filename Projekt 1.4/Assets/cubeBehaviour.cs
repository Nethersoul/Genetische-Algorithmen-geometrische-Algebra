using UnityEngine;
using System.Collections;

public class cubeBehaviour : MonoBehaviour {




	public Vector3 translation;
	public Vector3 rotation;
	public Vector3 scale = new Vector3(1, 1, 1);
	public float [,] matrix;

	// Use this for initialization
	void Start () {
		matrix = new float[,]{	{(float)Mathf.Cos(0.3f),	-((float)Mathf.Sin(0.3f)),	0.0f},
								{(float)Mathf.Sin(0.3f),	(float)Mathf.Cos(0.3f),		0.0f},
								{0.0f,						0.0f,						1.0f}
		};

	}
	
	// Update is called once per frame
	void Update () {


		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] vertices = mesh.vertices;
		Vector3[] normals = mesh.normals;
		int i = 0;

	


		while (i < vertices.Length) {


			vertices[i] += new Vector3(0, 10, 0);



			/*for(int j=0; j<3; j++){
				for (int s=0; s<3;s++){
					vertices[i]*=matrix[j,s];
//					Debug.Log("length is "+ matrix.GetLength(s)+ "and s is "+s);
				}
			}*/
//			vertices[i] += normals[i] * Random.value;
			i++;
		}
		mesh.vertices = vertices;
	
	}
}
