using System;

namespace GraphSynth.Search.Algorithms.Bandits {
    public class UCTBandit : AbstractBandit {

        private const double empiricalC = 1;

        public UCTBandit(int numArms) : base(numArms) {

        }

        /// <inheritdoc />
        /// <summary>
        /// Pulls the each arm once, then select an arm according to Pai(s) = argmax Q(s,a) + c * sqr(logn(s) / n(s,a))
        /// </summary>
        public override int SelectPullArm() {
            
            if (TotalPulls < NumArms) 
                return TotalPulls; // haven't tried each arm at least once yet

            var avgReward = GetAvgReward();
            var numPulls = GetNumPulls();
            var policyVal = new double[NumArms];

            double max = 0;
            var maxIndex = 0;
            for (var i = 0; i < NumArms; i++) {
                policyVal[i] = avgReward[i] + empiricalC * Math.Sqrt(Math.Log(TotalPulls) / numPulls[i]);
                if (policyVal[i] > max) {
                    max = policyVal[i];
                    maxIndex = i;
                }
            }
            return maxIndex;


        }
    }
}