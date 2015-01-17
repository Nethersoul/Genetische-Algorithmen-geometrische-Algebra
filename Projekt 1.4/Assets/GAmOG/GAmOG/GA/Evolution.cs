using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using GAmOG.utils;
using GAmOG.Tripod;

namespace GAmOG.GA
{
    public abstract class Evolution
    {
        protected Genom currentGenom = null;
        protected Genom futureGenom = null;
        private int mutationChance = 50;     //0-1000
        protected int NumbersOfCandidates = EvolutionControl.pop;

        private int trendSpares = 3;
        private int trendThreshold = 12;


        public Evolution()
        {
            this.defineGenom();
        }

        public virtual void evolve1()
        {
            for(int i=4;i<this.NumbersOfCandidates;i++)
            {
                DNA newChild = this.currentGenom.mate();
                newChild.mutate(mutationChance);
                futureGenom.addDNA(newChild);
            }
        }

		public virtual void evolve2()
		{
			this.futureGenom.calcFitnessTable();
			IList<int> fitness = this.futureGenom.getFitnessTable();
			this.calcMutationChance(fitness);
			
			//swap Genoms
			DNA best1 = this.futureGenom.getBestDNA ();
            DNA best2 = this.futureGenom.getDNA(1);
            DNA best3 = this.futureGenom.getDNA(2);
            DNA best4 = this.futureGenom.getDNA(3);

            UnityEngine.Debug.Log(best1.fitness );
            UnityEngine.Debug.Log(best2.fitness);
            UnityEngine.Debug.Log(best3.fitness);
            UnityEngine.Debug.Log(best4.fitness+" d");

			this.currentGenom = this.futureGenom;
			this.futureGenom = new Genom_Tripod(this.NumbersOfCandidates);

			this.futureGenom.setDNA(0,best1);
			this.futureGenom.setDNA(1, best2);
			this.futureGenom.setDNA(2, best3);
			this.futureGenom.setDNA(3, best4);

			//NOTE: futureGenom must be redefined in child class
		}

		public DNA_Walking getBestWalkDNA()
		{
			return (DNA_Walking)this.currentGenom.getBestDNA();
		}

		public DNA getDNA(int idx)
		{
			return this.futureGenom.getDNA (idx);
		}


		private void calcMutationChance(IList<int> fitness)
        {
			int trend = 0;
        //    int sumFitOld = 0;
            int sumFitNew = 0;


			IList<int> fitnessOld  = this.currentGenom.getFitnessTable();
            
            for(int i = 0; i<fitnessOld.Count;i++)
            {
                trend += (fitnessOld[i] - fitness[i]);
           //     sumFitOld += fitnessOld[i];
                sumFitNew += fitness[i];
            }



           if ((trend * trend) > this.trendThreshold * this.trendThreshold)
            {
                this.mutationChance -= 2;
                if(mutationChance<0)
                {
                    mutationChance =1;
                }
            }

            if ((trend * trend) < this.trendThreshold * this.trendThreshold)
            {
                this.trendSpares--;
            }

            if(trendSpares==0)
            {
                this.trendSpares = 3;
                this.mutationChance += 2;
            }

            if (sumFitNew < 50)
            {
           //     this.mutationChance += 2;
            }
            UnityEngine.Debug.Log("Mutation chance: " + this.mutationChance);
        }

        protected abstract void defineGenom();
    }
}
