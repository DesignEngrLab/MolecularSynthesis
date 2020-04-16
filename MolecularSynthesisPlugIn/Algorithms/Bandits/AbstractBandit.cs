using System;
using System.Linq;


namespace MolecularSynthesis.Algorithms 
{
    public abstract class AbstractBandit {
        protected readonly int NumArms;
        protected int TotalPulls;
        private readonly int[] _numPulls; // how many times we've pulled each arm
        private readonly double[] _averageReward;

        /// <summary>
        /// Initialize a new AbstractBandit.
        /// </summary>
        /// <param name="numArms">How many actions are available at the corresponding state.</param>
        protected AbstractBandit(int numArms) {
            NumArms = numArms;
            _averageReward = new double [NumArms]; // rewards initialized to 0
            _numPulls = new int [NumArms];
        }

        /// <summary>
        /// Update the given arm's pull count, its average reward, and the total pull count.
        /// </summary>
        /// <param name="arm">The index of the arm to update.</param>
        /// <param name="reward">The observed reward.</param>
        public void Update(int arm, double reward) {
            _numPulls[arm]++;
            _averageReward[arm] = (_averageReward[arm] * (_numPulls[arm] - 1) + reward) / _numPulls[arm];
            TotalPulls++;
        }

        /// <summary>
        /// Get the nunmer of pulls of each arm.
        /// </summary>
        public int[] GetNumPulls() {
            return _numPulls;  
        }
        
        
        /// <summary>
        /// Get the average reward
        /// </summary>
        public double[] GetAvgReward() {
            return _averageReward;  
        }
        
        
        
        /// <summary>
        /// Get the arm with the best average reward.
        /// </summary>
        public int GetBestArm() {
            return Array.IndexOf(_averageReward, GetBestReward());  // TODO makes two passes
        }

        /// <summary>
        /// Return the best average reward.
        /// </summary>
        private double GetBestReward() {
            return _averageReward.Max();
        }
        
        
        /// <summary>
        /// Print average reward and arm pulls of each option.
        /// </summary>
        public void PrintOptInfo() {
            for (var i = 0; i < _averageReward.Length; i++) {
                Console.WriteLine("Opt {0}: reward {1}, #tries {2}", i, _averageReward[i], _numPulls[i]);
            }
        }


        /// <summary>
        /// Return true if this Bandit has at least an arm to pull, else return false
        /// </summary>
        public bool HasOption() {
            return NumArms > 0;
        }


        /// <summary>
        /// Selects and returns a pull arm pursuant to the details of the bandit algorithm.
        /// </summary>
        public abstract int SelectPullArm();
        
        
    }
}