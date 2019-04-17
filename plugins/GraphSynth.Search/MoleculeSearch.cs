using System;
using System.IO;
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

            randomTrail();
            //mcts();
            //bfs();

        }

        private void bfs()
        {
            var agent = new BFS(settings);
            Seed = new candidate(OBFunctions.tagconvexhullpoints(settings.seed), settings.numOfRuleSets);
            agent.search(Seed.copy(), 5);
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

        private void randomTrail() {
            var agent = new Algorithms.Random(settings);
            Seed = new candidate(OBFunctions.tagconvexhullpoints(settings.seed), settings.numOfRuleSets);
            agent.buildTrail(Seed, 5, 10);
        }


        
        

        public override string text => "Molecule Search Runner";
    }
}