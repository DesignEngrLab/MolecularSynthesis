using System;
using GraphSynth.Representation;
using GraphSynth.Search.Algorithms;
using OpenBabelFunctions;



namespace GraphSynth.Search {
    /// <summary>
    /// Enables plug-and-play comparison of algorithms.
    /// </summary>
    public class MoleculeSearch : SearchProcess {
        private candidate Seed;


        /// <inheritdoc />
        /// <summary>
        /// Initializes SearchProcess properties.
        /// </summary>
        public MoleculeSearch() {
            RequireSeed = true;
            RequiredNumRuleSets = 1;
            AutoPlay = true;

        }

        protected override void Run() {

            //mcts();
            //bfs();

        }



        private void mcts() {
            var time0 = DateTime.Now;
            var agent = new Mcts(settings);
            agent.loadSmileAndReward();
            Seed = new candidate(OBFunctions.tagconvexhullpoints(settings.seed), settings.numOfRuleSets);
            agent.buildTree(Seed.copy());
            agent.saveSmileAndReward();
            var time1 = DateTime.Now;
            Console.WriteLine("Time used: {0}", (time1 - time0).TotalSeconds);
            
            
            
        }
        public override string text => "Molecule Search Runner";
    }
}