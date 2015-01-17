using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GAmOG.GA
{
   public abstract class Gen
   {

       public abstract Gen crossover(Gen g, int crossoverPoint);

       public abstract void mutate();
   }
}
