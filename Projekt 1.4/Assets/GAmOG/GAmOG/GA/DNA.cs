using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GAmOG.GA
{
    public abstract class DNA :IComparable
    {
        protected List<Gen> Genes;
		public int fitness { get;protected set; }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="chance">0-100 The % chance to mutate</param>
        public abstract void  mutate(int chance);

        public abstract void calcFitness();


        public abstract DNA crossover(DNA parent2);

        public Gen getGene(int idx)
        {
            return this.Genes[idx];
        }



        public void addGen(Gen g)
        {
  
            this.Genes.Add(g);
        }

        public void clearGenes()
        {

            this.Genes.Clear();
        }

        public int CompareTo(DNA obj)
        {
            if(this.fitness==obj.fitness)
            {
                return 0;
            }
            else if(this.fitness>obj.fitness)
            {
                return -1;
            }
            else
            {
                return 1;
            }
        }

        public int CompareTo(object obj)
        {
            try
            {
                return this.CompareTo((DNA)obj);
            }
            catch(Exception x)
            {
               Console.WriteLine(" Error not implemented: Compare obj with DNA in DNA compare to msg: "+x.Message);
            }
            return 0;
        }
    }


}
