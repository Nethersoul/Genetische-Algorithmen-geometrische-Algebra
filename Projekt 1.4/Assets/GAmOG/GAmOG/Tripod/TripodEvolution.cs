using GAmOG.GA;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GAmOG.Tripod
{
    public class TripodEvolution : Evolution
    {

        public TripodEvolution():base()
        {
        
        }

        protected override void defineGenom()
        {
            this.currentGenom = new Genom_Tripod(this.NumbersOfCandidates);
            this.futureGenom  = new Genom_Tripod(this.NumbersOfCandidates);
        }

		public int getBestFitnessScore()
        {
            return this.currentGenom.getFitnessTable()[0];
        }

        public override void evolve1()
        {
            base.evolve1();
        }

		public override void evolve2()
		{

			base.evolve2();
		

		//	this.futureGenom = new Genom_Tripod(this.NumbersOfCandidates);
		
		}
    }
}
