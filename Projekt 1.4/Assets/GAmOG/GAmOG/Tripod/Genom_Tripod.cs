using GAmOG.GA;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using GAmOG.utils;

namespace GAmOG.Tripod
{
    /// <summary>
    /// This class holds all possible DNA of a Tribod. 
    /// </summary>
    public class Genom_Tripod : Genom
    {

        public Genom_Tripod(int amountCandidates)
        {
            this.DNACandidates = new List<DNA>(amountCandidates);
            this.defineDNA(amountCandidates);
        }

        private void defineDNA(int amountCandidates)
        {
            for(int i = 0; i<amountCandidates;i++)
            {
                DNA_Walking c = new DNA_Walking();
                this.DNACandidates.Add(c);
            }

        }

        /// <summary>
        /// Works only if fitness has been calculated
        /// </summary>
        /// <param name="parent1"></param>
        /// <param name="parent2"></param>
        /// <returns></returns>
        public override DNA crossover(DNA parent1, DNA parent2)
        {
         
            DNA newChild = parent1.crossover(parent2);
            return newChild;
        }

      
    }
}
