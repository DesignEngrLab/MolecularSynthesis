using System.Collections.Generic;
using GraphSynth;
using GraphSynth.Representation;
using GraphSynth.Search;
using Priority_Queue;
using PropertyEvaluation;
namespace GearSynth
{
    public class GearSearch : SearchProcess
    {
        // give desiredMoment
        double[] desiredMoment = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

        public GearSearch(GlobalSettings settings) : base(settings)
        {
        }

        /// <summary>
        ///   Gets the text describing that is displayed in the menu. It must be overridden in the methods.
        /// </summary>
        /// <value>The text.</value>
        public override string text
        {
            get { return "Run Gear Synth"; }
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
                var options = rulesets[0].recognize(current.graph);


                for (var i = 0; i != options.Count; i++)
                {
                    var child = current.copy();
                    transferLmappingToChild(child.graph, current.graph, options[i]);
                    options[i].apply(child.graph, null);
                    child.addToRecipe(options[i]);
                    //SearchIO.output("...done.", 3);
                    //need to have the desird moment and the distance function
                    //child.f3 = Evaluation.CalcMoment(child);
                    child.f0 = Evaluation.distance(child, desiredMoment);
                    candidates.Enqueue(child, child.f0);
                    System.Console.WriteLine(child.f0);


                    /* f0 is mass/weight in lbs.
                     * f1 is inefficiency. 0 is 100% efficient, and 100 is 95 is 5% efficient
                     * f2 is the amount of constraint violation.
                     * f3 is a weighted sum of these two based on value in UserDefinedGoals */
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
            /* f2 represents the amount of infeasibility in the design. We are seeking
             * the first design that is feasible. Following the negative-null form
             * convention, feasibility occurs at zero or below. */
            return (current.f2 <= 0);
        }
    }
}