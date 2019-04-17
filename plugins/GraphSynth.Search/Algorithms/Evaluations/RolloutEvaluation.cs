using System;
using System.IO;
using System.Collections.Generic;
using GraphSynth.Representation;
using OpenBabelFunctions;

namespace GraphSynth.Search.Algorithms.Evaluations {
    public class RolloutEvaluation : AbstractEvaluation {
        private Random _random;
        
        public RolloutEvaluation(int width, int depth, AbstractAlgorithm parent) : base(width, depth, parent) {
            _random = new Random(Parent.Settings);
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
            
            var carboxylOpts = Parent.GetCarboxylOptions(candCopy);
            if (carboxylOpts.Count != 0) {
                for (var i = 0; i < carboxylOpts.Count; i++) {
                    var evaluateCopy = candCopy.copy();
                    var evaluateOpts = Parent.GetCarboxylOptions(evaluateCopy);
                    Parent.ApplyOption(evaluateOpts[i], evaluateCopy, true);
                    var reward = Parent.Evaluate(evaluateCopy);
                    if (reward > bestRewardSoFar) {
                        bestRewardSoFar = reward;
                    }
                }

            }
            
            for (var iteration = 0; iteration < Depth; iteration++) {
                
                var opt = _random.ChooseOption(candCopy);
                if (opt == null)
                    break;
                Parent.ApplyOption(opt, candCopy, true);
                
                carboxylOpts = Parent.GetCarboxylOptions(candCopy);
                if (carboxylOpts.Count != 0) {
                    for (var i = 0; i < carboxylOpts.Count; i++) {
                        var evaluateCopy = candCopy.copy();
                        var evaluateOpts = Parent.GetCarboxylOptions(evaluateCopy);
                        Parent.ApplyOption(evaluateOpts[i], evaluateCopy, true);
                        var reward = Parent.Evaluate(evaluateCopy);
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