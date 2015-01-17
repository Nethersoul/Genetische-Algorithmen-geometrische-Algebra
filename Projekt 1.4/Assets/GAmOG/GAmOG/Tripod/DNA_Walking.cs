using GAmOG.GA;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using GAmOG.utils;

namespace GAmOG.Tripod
{
    public class DNA_Walking : DNA
    {
        
        public DNA_Walking()
        {
            base.Genes = new List<Gen>(EvolutionControl.pop);
            this.defineGenes();
        }

        private void defineGenes()
        {
            int JointsAngleFrom_a3d = 80;       //80
            int JointsAngleFrom_b3d = 0;        //0
            int JointsAngleTo_a3d = 120;         //120
            int JointsAngleTo_b3d = 30;         //30

            int JointsAngleFrom = 30;        //0
            int JointsAngleTo = 120;        //120

            //Artificial splitting of 3d joints in 2 parts
            Gen_Joint joint1_a = new Gen_Joint(0, 1, 0, JointsAngleFrom_a3d, JointsAngleTo_a3d);
            Gen_Joint joint1_b = new Gen_Joint(0, 0, -1, JointsAngleFrom_b3d, JointsAngleTo_b3d);		//0.0.1


            Gen_Joint joint2_a = new Gen_Joint(0, -0.5, -0.5, JointsAngleFrom_a3d, JointsAngleTo_a3d);
            Gen_Joint joint2_b = new Gen_Joint(0, -0.5, 0.5, JointsAngleFrom_b3d, JointsAngleTo_b3d);       //0/-0.5/0.5

            Gen_Joint joint3_a = new Gen_Joint(0, -0.5, 0.5, JointsAngleFrom_a3d, JointsAngleTo_a3d);
            Gen_Joint joint3_b = new Gen_Joint(0, 0.5,0.5, JointsAngleFrom_b3d, JointsAngleTo_b3d);     

            //3 lower joint
            Gen_Joint joint4 = new Gen_Joint(0, 1, 0, JointsAngleFrom, JointsAngleTo);
            Gen_Joint joint5 = new Gen_Joint(0, -0.5, -0.5, JointsAngleFrom, JointsAngleTo);
            Gen_Joint joint6 = new Gen_Joint(0, -0.5, 0.5, JointsAngleFrom, JointsAngleTo);

            //one lowest joint
            Gen_Joint joint7 = new Gen_Joint(0, 1, 0, JointsAngleFrom, JointsAngleTo);

        
            if (this.Genes.Count == 0)
            {
                this.Genes.Add(joint1_a);
                this.Genes.Add(joint1_b);
                this.Genes.Add(joint2_a);
                this.Genes.Add(joint2_b);
                this.Genes.Add(joint3_a);
                this.Genes.Add(joint3_b);
                this.Genes.Add(joint4);
                this.Genes.Add(joint5);
                this.Genes.Add(joint6);
                this.Genes.Add(joint7);
            }
        }

        public override void mutate(int chance)
        {
            foreach(Gen g in this.Genes)
            {
                if (Rand.Instance.randomNumber(1000) < (chance))
                {
                    g.mutate();
                }
            }
        }

        override
        public void calcFitness()
        {
           //TODO: this
           // this.fitness = 42;
        }

		public void defineFitness(int fit)
		{
			this.fitness = fit;
		}


        public override DNA crossover(DNA parent2)
        {
            int crossoverPoint = Rand.Instance.randomNumber(RotorTimeMutant.nrOfChanges); //TODO replace 9
            DNA_Walking newChild = new DNA_Walking();
            newChild.clearGenes();

            for(int i = 0; i<this.Genes.Count; i++)
            {   
                newChild.addGen(this.Genes[i].crossover(parent2.getGene(i),crossoverPoint));
            }

            return newChild;
        }
    }
}
