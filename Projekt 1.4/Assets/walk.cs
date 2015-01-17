using UnityEngine;
using System.Collections;
using GAmOG.Tripod;
using c3ga_ns;

public class walk : MonoBehaviour 
{

	public GameObject Leg1;
	public GameObject Leg2;
	public GameObject Leg3;
	public GameObject UpperLeg1;
	public GameObject UpperLeg2;
	public GameObject UpperLeg3;
	public GameObject UpperUpperLeg1;
	

	public DNA dna;
	private DNA_Walking d;

	public struct DNA
	{
		public Gen_Joint Leg1DNA;
		public Gen_Joint Leg2DNA;
		public Gen_Joint Leg3DNA;
		public Gen_Joint UpperLeg1DNA;
		public Gen_Joint UpperLeg2DNA_a;
		public Gen_Joint UpperLeg2DNA_b;
		public Gen_Joint UpperLeg3DNA_a;
		public Gen_Joint UpperLeg3DNA_b;
		public Gen_Joint UpperUpperLeg1DNA_a;
		public Gen_Joint UpperUpperLeg1DNA_b;

		public DNA(DNA_Walking walkingDna)
		{
			this.Leg1DNA = (Gen_Joint)walkingDna.getGene (9);
			this.UpperLeg1DNA = (Gen_Joint)walkingDna.getGene (6);
			this.UpperUpperLeg1DNA_a = (Gen_Joint)walkingDna.getGene (0);
			this.UpperUpperLeg1DNA_b = (Gen_Joint)walkingDna.getGene (1);
			
			this.Leg2DNA = (Gen_Joint)walkingDna.getGene (7);
			this.UpperLeg2DNA_a = (Gen_Joint)walkingDna.getGene (2);
			this.UpperLeg2DNA_b = (Gen_Joint)walkingDna.getGene (3);

			this.Leg3DNA = (Gen_Joint)walkingDna.getGene (8);
			this.UpperLeg3DNA_a = (Gen_Joint)walkingDna.getGene (4);
			this.UpperLeg3DNA_b = (Gen_Joint)walkingDna.getGene (5);

		}
	}




	// Use this for initialization
	void Start ()
	{

	
	}

	private Quaternion readQuat(Gen_Joint d)
	{
		double[] quats = d.getQuaternion(Mathf.PingPong (Time.time*10, 5));
		return new Quaternion ((float)quats [0], (float)quats [1], (float)quats [2], (float)quats [3]);
	}

	private void transformLeg(GameObject Leg, Gen_Joint dna)
	{
		Leg.transform.localRotation = readQuat (dna);//new Quaternion ();// readQuat (dna);
	}

	private void transformLeg(GameObject Leg, Gen_Joint dna_a, Gen_Joint dna_b, Quaternion minusOffset_b)
	{
	//	Debug.Log (readQuat (dna_a).eulerAngles);
		Leg.transform.localRotation = minusOffset_b*readQuat (dna_b)*readQuat (dna_a);//Quaternion.Euler(-90,0,0)*Quaternion.Euler(0,90,0)*readQuat (dna_a); //new Quaternion ();
	}

	public void defineFit(int fit)
	{
		this.d.defineFitness (fit);
	}

	public void setWalkdingDna(DNA_Walking d )
	{
		this.d = d;
	}


	void FixedUpdate() 	//FixedU
	{
		if (this.d != null) {
						transformLeg (Leg1, dna.Leg1DNA);
						transformLeg (UpperLeg1, dna.UpperLeg1DNA);
						transformLeg (UpperUpperLeg1, dna.UpperUpperLeg1DNA_a, dna.UpperUpperLeg1DNA_b, Quaternion.Euler (0, 0, 15));
			
						transformLeg (Leg2, dna.Leg2DNA);	
						transformLeg (UpperLeg2, dna.UpperLeg2DNA_a, dna.UpperLeg2DNA_b, Quaternion.Euler (0, 7.5f, -7.5f));

						transformLeg (Leg3, dna.Leg3DNA);
						transformLeg (UpperLeg3, dna.UpperLeg3DNA_a, dna.UpperLeg3DNA_b, Quaternion.Euler (0, -7.5f, -7.5f));
				} 
		
	}


	// Update is called once per frame

}
