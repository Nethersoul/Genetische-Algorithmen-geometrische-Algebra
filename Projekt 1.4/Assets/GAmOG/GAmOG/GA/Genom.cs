using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using GAmOG.utils;

namespace GAmOG.GA
{
    public abstract class Genom
    {
        protected List<DNA> DNACandidates;
 
		protected int SumFitness = 0;


        public abstract DNA crossover(DNA parent1, DNA parent2);


        /// <summary>
        /// Selection of parents based on fitness
        /// </summary>
        /// <returns></returns>
        public DNA mate()
        {
            int NDna = this.DNACandidates.Count;

            int idx1 = 0;
            int idx2 = 0;
            int i = 0;
            int selectionNumber = Rand.Instance.randomNumber(SumFitness);
            var keys = from dna in this.DNACandidates select dna.fitness;

           

            for (i = 0; selectionNumber >= 0; i++)
            {
                if (i ==( keys.Count()-1))
                {
                    break;
                }
                int k = keys.ElementAt<int>(i);
                selectionNumber -= k;
            }
            idx1 = i;

            //UnityEngine.Debug.Log(selectionNumber + " " + SumFitness);
         
           selectionNumber = Rand.Instance.randomNumber(SumFitness);
            for (i = 0; selectionNumber >= 0; i++)
            {
                if (i == (keys.Count() - 1))
                {
                    break;
                }
                int k = keys.ElementAt<int>(i);
                selectionNumber -= k;
            }

            idx2 = i;

            return this.crossover(this.DNACandidates[idx1], this.DNACandidates[idx2]);
        }

        public void calcFitnessTable()
        {
            SumFitness = 0;
           
            foreach (DNA d in this.DNACandidates)
            {
                d.calcFitness();
                SumFitness += d.fitness;
               
            }
            this.DNACandidates.Sort();
        }

        /// <summary>
        /// Call calcFitnessTable before this method
        /// </summary>
        /// <returns>Null if no fitness is know, Best DNA otherwise</returns>
        public  DNA getBestDNA()
        {
            var x = from d in this.DNACandidates where d.fitness == this.DNACandidates.First().fitness select d;
           

            int idx = Rand.Instance.randomNumber(x.Count() - 1);
            foreach (DNA a in x)
            {
                if(idx==0)
                {
                    return a;
                }
                idx--;
            }

            return null;
        }

		public void setDNA(int i, DNA d)
		{
			this.DNACandidates [i] = d;
		}

	


		public DNA getDNA(int i)
		{
			return this.DNACandidates [i];
		}


        public void addDNA(DNA a)
        {
            this.DNACandidates.Add(a);
        }

		public IList<int> getFitnessTable()
        {
            var x = from d in this.DNACandidates select d.fitness;
			foreach(DNA d in this.DNACandidates)

			{
				int a = d.fitness;
				d.calcFitness();
			}
            return x.ToList();
        }
    }
}
