using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace GAmOG.utils
{
    public class Rand
    {
        private static Rand instance;
        private static Random r;

        private Rand() { }

        public static Rand Instance
        {
            get
            {
                if (instance == null)
                {
                    instance = new Rand();
                    r = new Random();
                }
                return instance;
            }
        }

        public int randomNumber(int maxRange)
        {
            return r.Next(maxRange);
        }
    }

}
