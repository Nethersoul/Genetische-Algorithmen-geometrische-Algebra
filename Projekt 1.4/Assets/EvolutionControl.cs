using UnityEngine;
using System.Collections;
using System.Collections.Generic;
using GAmOG.Tripod;


public class EvolutionControl : MonoBehaviour
{
	
	private TripodEvolution evo = new TripodEvolution();
	ArrayList tripods = new ArrayList();

	private double walktimeDef = 30;
	private double walktime ;

	private bool ready = false;
    public static int pop = 15;
	
	void evolve1()
	{
		evo.evolve1();
        for (int i = 0; i < pop; i++)
		{
		
			GameObject triCopy = (GameObject)this.tripods[i];
			triCopy.GetComponent<walk>().dna = new walk.DNA((DNA_Walking)this.evo.getDNA(i));
		
			triCopy.GetComponent<walk>().setWalkdingDna((DNA_Walking)this.evo.getDNA(i));

		}

	}

	void evolve2()
	{
        for (int i = 0; i < pop; i++) 
		{
			GameObject tri = (GameObject)this.tripods[i];
			int fit =(int)tri.GetComponent<Transform>().position.z;
            if (fit < 0)
            {
                fit = 0;
            }
			tri.GetComponent<walk>().defineFit(fit);
		}



		evo.evolve2();

		Debug.Log ("BestScore" + evo.getBestFitnessScore () + "\n");

	}

	// Use this for initialization
	void Start () 
	{
		walktime = walktimeDef;
		GameObject triOrg  = GameObject.Find ("Tripod");

        for (int i = 0; i < pop; i++)
		{
			GameObject triCopy = (GameObject)Instantiate(triOrg);
			triCopy.transform.position = new Vector3(i*60-300,0,0);
			triCopy.transform.rotation = new Quaternion();
			triCopy.GetComponent<walk>().dna =  new walk.DNA((DNA_Walking)this.evo.getDNA(i));
			triCopy.GetComponent<walk>().setWalkdingDna((DNA_Walking)this.evo.getDNA(i));
   			this.tripods.Add(triCopy);
            
		}
		triOrg.transform.position = new Vector3 (0, 0, -500);

	//	Physics.gravity = new Vector3 (0, -100, 0);
		ready = true;
	}
	
	// Update is called once per frame
	void Update () 
	{
		if (ready)
		{
			walktime -= Time.deltaTime;
			if (walktime <= 0) 
			{
				evolve2 ();
				evolve1 ();
                for (int i = 0; i < pop; i++)
				{
					GameObject tri = (GameObject)this.tripods[i];
					tri.transform.position = new Vector3(i*60-300,3,0);
					tri.transform.rotation = new Quaternion();
				}
				this.walktime = walktimeDef;
			}
		}
	}
}
