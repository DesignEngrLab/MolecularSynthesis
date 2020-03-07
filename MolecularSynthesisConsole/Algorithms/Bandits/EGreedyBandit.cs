using System;

namespace GraphSynth.Search.Algorithms.Bandits {
    public class EGreedyBandit : AbstractBandit {
        private readonly double _epsilon;
        private readonly System.Random _rand = new System.Random();

        public EGreedyBandit(int numArms, double epsilon) : base(numArms) {
            _epsilon = epsilon;
        }

        /// <inheritdoc />
        /// <summary>
        /// Pulls the most rewarding arm with (1 - epsilon) probability; otherwise, a non-optimal arm is pulled at random.
        /// </summary>
        public override int SelectPullArm() {
            
            if (TotalPulls < NumArms) return TotalPulls; // haven't tried each arm at least once yet

            var bestArm = GetBestArm();
            if (!(_rand.NextDouble() < _epsilon)) return bestArm;
            
            var nonBest = _rand.Next(NumArms);
            while (nonBest == bestArm)
                nonBest = _rand.Next(NumArms);
            return nonBest; // pull a random non-best arm
        }
    }
}