using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using c3ga_ns;
using GAmOG.utils;

namespace GAmOG.Tripod
{
    /// <summary>
    /// This class is responsible for changes of the angle. The rotor just defines where and how much that can be turned. 
    /// This class does define by how much that shall be turned at a certain timepoint.
    /// </summary>
    public class RotorTimeMutant
    {
        public static int nrOfChanges = 18;       //at max 10 changes are done in one iteration
        private float[] changefactors;    //0 to 1 turns
        public RotorTimeMutant()
        {
            this.changefactors = new float[RotorTimeMutant.nrOfChanges];
            for (int i = 0; i < RotorTimeMutant.nrOfChanges; i++)
            {
                //initialisation with 0.1 to 1 meaning one full turn over the joint
                this.changefactors[i] = (i+1)/(EvolutionControl.pop*1f);//((Rand.randomNumber(2000) - 1000) / 1000f);
            }
        }

        public RotorTimeMutant(RotorTimeMutant rtm)
        {
            this.changefactors = new float[RotorTimeMutant.nrOfChanges];
            for (int i = 0; i < RotorTimeMutant.nrOfChanges; i++)
            {
                this.changefactors[i] = rtm.changefactors[i];
            }
        }

        /// <summary>
        /// Interpolates angle depending on which timeframe. If we have a cycle of 100 frames and are on the 6 frame 
        /// than the value 6/100*10= 0.6 is given to this function (*10 since from 0-10) where 10 = nrOfChanges
        /// 
        /// </summary>
        /// <param name="time">Time between 0-(nrOfChanges-1) in float</param>
        /// <returns>Factor which will be multiplied with angle</returns>
        public float getChangeFactor(float time)
        {
            int i = (int)Math.Floor(time);
            if(i==(RotorTimeMutant.nrOfChanges-1)||(time-i)<0.01)
            {
                return this.changefactors[i];
            }

            float w1 = this.changefactors[i];
            float w2 = this.changefactors[i + 1];

            return w1 + (w2 - w1) *  (time-i);
        }

        public void getCrossed(RotorTimeMutant rm, int point)
        {
            for (int i = point; i < RotorTimeMutant.nrOfChanges; i++)
            {
                this.changefactors[i] = rm.changefactors[i];
            }
        }

        public void getMutated()
        {
            for (int i = 0; i < RotorTimeMutant.nrOfChanges; i++)
            {
                //chance for a mutation is about 10% for each-> might not mutate at all or do multiple mutation at once
                if (Rand.Instance.randomNumber(RotorTimeMutant.nrOfChanges - 1) == 1)
                {
                    this.changefactors[i] = (Rand.Instance.randomNumber(1000)) / 1000f;   //number in float from 0 to +1
                }
            }
        }

    }
}
