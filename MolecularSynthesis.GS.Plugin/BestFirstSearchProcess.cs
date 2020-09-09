using System.Collections.Generic;
using System.Threading;
using GraphSynth;
using GraphSynth.Representation;
using GraphSynth.Search;
using Priority_Queue;

namespace MolecularSynthesis.Plugin
{
    public class BestFirstSearch : SearchProcess
    {
        // give desiredMoment
        // [] desiredMoment = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        double[] desiredLenghtAndRadius = new double[] { 300, 50 };

        public BestFirstSearch(GlobalSettings settings) : base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 2;
            AutoPlay = true;
        }

        /// <summary>
        ///   Gets the text describing that is displayed in the menu. It must be overridden in the methods.
        /// </summary>
        /// <value>The text.</value>
        public override string text
        {
            get { return "Best First Search"; }
        }

        /// <summary>
        ///   Runs the instance of the plugin that is initialized by GraphSynth.
        ///   This is basically a best-first tree search for the first solution
        ///   that doesn't violate constraints.
        /// </summary>
        protected override void Run()
        {
            var candidates = new SimplePriorityQueue<candidate, double>();
            var current = seedCandidate;
            candidates.Enqueue(current, 0);
            //candidates.Enqueue(current, -1);
            int iteration = 0;
            while (!SearchIO.terminateRequest && (candidates.Count != 0))
            {
                current = candidates.Dequeue();
                if (isCurrentTheGoal(current)) break;

                //SearchIO.output("Current: mass=" + current.f0
                                //+ " ineff=" + current.f1 + " constraintViolation=" + current.f2);
                SearchIO.iteration = iteration++;
                SearchIO.miscObject = current.recipe.Count;
                //RECOGNIZE
                //var options = rulesets[0].recognize(current.graph);
                //rulesets[0].nextRuleSet(GenerationStatuses.Unspecified);
                //current.activeRuleSetIndex;
                
                var childrenCandidate=RecognizeChooseApply.GenerateAllNeighbors(current,rulesets,false,false,true);
                foreach (var child in childrenCandidate)
                {
                    //child.f0 = Evaluation.distance(child, desiredLenghtAndRadius);
                    candidates.Enqueue(child, child.f0);
                    SearchIO.output(child.f0,3);
                }
                
            }
            SearchIO.addAndShowGraphWindow(current.graph, "Here is your gear train.");
            Save(settings.OutputDir + "testTuned", current);
        }


        /// <summary>
        ///   Determines whether [is current the goal] [the specified current].
        /// </summary>
        /// <param name = "current">The current.</param>
        /// <returns>
        ///   <c>true</c> if [is current the goal] [the specified current]; otherwise, <c>false</c>.
        /// </returns>
        private static bool isCurrentTheGoal(candidate current)
        {
         
            return (current.f2 <= 0);
        }
    }
}