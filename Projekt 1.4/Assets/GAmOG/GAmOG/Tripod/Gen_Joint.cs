using GAmOG.GA;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using c3ga_ns;
namespace GAmOG.Tripod
{
    public class Gen_Joint : Gen
    {
        double from = 0;
        private rotor r;
        private RotorTimeMutant rm;


        /// <summary>
        /// Erstellt einen Rotor mit dem Skalar 0 und der Ebene definiert durch e1^e2, e2^e3 und e4^e1
        /// </summary>
        /// <param name="e1_e2">x-y ebene</param>
        /// <param name="e2_e3">y-z ebene</param>
        /// <param name="e3_e1">z-x ebene</param>
        public Gen_Joint(double e1_e2, double e2_e3, double e3_e1, float angleFrom, float angleTo)
        {
            double radFrom = (Math.PI / 180) * angleFrom;
            double radTo = (Math.PI / 180) * angleTo ;
           
            double phi = (utils.Rand.Instance.randomNumber((int)((radTo*100 - radFrom*100))) / 100f);      //amount of turn given in rad (means from 30° - 180° is possible)
            this.from = radFrom;
         
            
            // UnityEngine.Debug.Log(radFrom + " " + radTo);

            double cFac = Math.Cos(phi*2);

            r = new rotor(rotor.CoordinateOrder.coord_scalar_e1e2_e2e3_e3e1, Math.Sin(phi*2), e1_e2 * cFac, e2_e3 * cFac, e3_e1 * cFac);
            this.rm = new RotorTimeMutant();
        }

        public Gen_Joint(rotor r,RotorTimeMutant rm)
        {
            this.r = new rotor(r);
            this.rm = new RotorTimeMutant(rm);
        }

        public override Gen crossover(Gen g, int crossoverPoint)
        {
           return this.crossover((Gen_Joint)g, crossoverPoint);
        }

        public Gen crossover(Gen_Joint g, int crossoverPoint)
        {
            Gen_Joint gNew = new Gen_Joint(this.r,this.rm);
            gNew.crossover(g.getRotorTimeMutant(),crossoverPoint);
            return gNew;
        }


        public void crossover(RotorTimeMutant rtm,int  crossoverPoint)
        {
            this.rm.getCrossed(rtm, crossoverPoint);
        }

        public RotorTimeMutant getRotorTimeMutant()
        {
            return this.rm;
        }

        /// <summary>
        /// Rotates this rotor with a given rotor r
        /// </summary>
        /// <param name="r">Rotor to apply on this rotor</param>
        public void rotate(rotor r)
        {
            this.r = c3ga.gp(this.r, r);
        }

        public double[] getQuaternion(float time)
        {
			float f = this.getRotorTimeMutant ().getChangeFactor (time);

            double sinPhi = this.r.get_scalar();
			double phi = Math.Asin(sinPhi)*2*f+this.from;
       
         //   phi = this.from;
			double cosPhi = Math.Cos(phi/2);        //cause in rotor phi is doubled
			sinPhi = Math.Sin(phi/2);

            double[] rots = new double[4];
        

            rots[0] = this.r.get_e1_e2() / cosPhi * sinPhi;
            rots[1] = this.r.get_e2_e3() / cosPhi * sinPhi;
            rots[2] = this.r.get_e3_e1() / cosPhi * sinPhi;
            rots[3] = cosPhi;
            return rots;
        }


        public override void mutate()
        {
            this.rm.getMutated();
        }
    }
}
