using System;
using System.IO;
using System.Collections.Generic;
using GraphSynth.Representation;
using OpenBabelFunctions;

namespace GraphSynth.Search.Algorithms.Evaluations {
    public class RolloutEvaluation : AbstractEvaluation {
        private Random _random;
        
        public RolloutEvaluation(int width, int depth, AbstractAlgorithm parent) : base(width, depth, parent) {
            _random = new Random(AbstractAlgorithm.Settings);
        }


        public override double Evaluate(candidate state) {
            double reward = 0;
            
            for (var i = 0; i < Width; i++)  // do Width simulations
                reward += RunRolloutLocal(state);
            return reward / Width;
        }


        /// <summary>
        /// Do one rollout for the state.
        /// </summary>
        private double RunRolloutLocal(candidate cand) {
            var candCopy = cand.copy();
            double bestRewardSoFar = 0;
            
            var carboxylOpts = AbstractAlgorithm.GetCarboxylOptions(candCopy);
            if (carboxylOpts.Count != 0) {
                for (var i = 0; i < carboxylOpts.Count; i++) {
                    var evaluateCopy = candCopy.copy();
                    var evaluateOpts = AbstractAlgorithm.GetCarboxylOptions(evaluateCopy);
                    AbstractAlgorithm.ApplyOption(evaluateOpts[i], evaluateCopy, true);
                    var reward = AbstractAlgorithm.Evaluate(evaluateCopy);
                    if (reward > bestRewardSoFar) {
                        bestRewardSoFar = reward;
                    }
                }

            }
            
            for (var iteration = 0; iteration < Depth; iteration++) {
                
                var opt = _random.ChooseOption(candCopy);
                if (opt == null)
                    break;
                AbstractAlgorithm.ApplyOption(opt, candCopy, true);
                
                carboxylOpts = AbstractAlgorithm.GetCarboxylOptions(candCopy);
                if (carboxylOpts.Count != 0) {
                    for (var i = 0; i < carboxylOpts.Count; i++) {
                        var evaluateCopy = candCopy.copy();
                        var evaluateOpts = AbstractAlgorithm.GetCarboxylOptions(evaluateCopy);
                        AbstractAlgorithm.ApplyOption(evaluateOpts[i], evaluateCopy, true);
                        var reward = AbstractAlgorithm.Evaluate(evaluateCopy);
                        if (reward > bestRewardSoFar) {
                            bestRewardSoFar = reward;
                        }
                    }

                }


            }
            return bestRewardSoFar;
        }
    }
}