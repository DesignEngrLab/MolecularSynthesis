using GraphSynth.Representation;

namespace MolecularSynthesis
{
    /// <summary>
    /// The main class for implementing state evaluation functions.
    /// </summary>
    public abstract class AbstractEvaluation {
        protected int Width;
        protected int Depth;
        protected AbstractAlgorithm Parent;  // how we can call f
        
        /// <summary>
        /// Implements a function for performing random rollout to the given parameters.
        /// </summary>
        /// <param name="width">How many rollouts to run.</param>
        /// <param name="depth">How deep each rollout should go.</param>
        /// <param name="parent">The class used to manipulate the molecule.</param>
        protected AbstractEvaluation(int width, int depth, AbstractAlgorithm parent) {
            Width = width;
            Depth = depth;
            Parent = parent;
        }

        /// <summary>
        /// Using the predefined policy/policies, width, and depth, evaluates the state.
        /// </summary>
        /// <param name="state"></param>
        /// <returns></returns>
        public abstract double Evaluate(candidate state);
    }
}