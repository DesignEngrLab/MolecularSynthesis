using GraphSynth.Representation;
using LAMMPSnow;
using OpenBabel;
using OpenBabelFunctions;
using System;
using System.Collections.Generic;
using System.IO;
using library;



namespace GraphSynth.Search.Algorithms {
    
    public abstract class AbstractAlgorithm {
        public readonly GlobalSettings Settings;
        public readonly System.Random Rand = new System.Random();
        private static readonly string IODir = OBFunctions.GetRamDir();

        public readonly string _runDirectory;
        public abstract string RunDirectoryName { get; } // the algorithm-based output directory
        private readonly string _inputFilePath;

        private const double Slope = 0.5; // kcal/mol/angstrom over strain percent
        private const double AngleFloor = 155; // minimum acceptable angle between carboxylates

        public static readonly graph2almostanything Converter = new graph2almostanything();
        private readonly Dictionary<string, double> previous = new Dictionary<string, double>(); // cache evaluations
        private readonly HashSet<string> uniqueLinkers = new HashSet<string>();
        private int candidateCnt = 0;
        


        protected AbstractAlgorithm(GlobalSettings settings_) {
            Settings = settings_;

            // Make sure we can write output files
            _runDirectory = Path.Combine(Settings.OutputDirAbs, RunDirectoryName);
            //_inputFilePath = Path.Combine(Settings.OutputDirAbs, "config/in.deform");
            if (!Directory.Exists(_runDirectory))
                Directory.CreateDirectory(_runDirectory);
            //if (!File.Exists(_inputFilePath))
            //    throw new Exception("Input file doesn't exist: " + _inputFilePath);
        }
        
        public void loadSmileAndReward (){
            if (File.Exists(_runDirectory + "/dictInfo.txt")) {
                using (StreamReader file = new StreamReader(_runDirectory + "/dictInfo.txt")) {
                    string line;
                    char[] separator = {','};
                    while ((line = file.ReadLine()) != null) {
//                    Console.WriteLine(line);
                        string[] splits = line.Split(separator, StringSplitOptions.None);
                        previous.Add(splits[0], Convert.ToDouble(splits[1]));
                        candidateCnt++;
                    }
                }
            }
        }

        public void saveSmileAndReward (){
            using (StreamWriter file = new StreamWriter(_runDirectory + "/dictInfo.txt")) {
                foreach (var entry in previous) {
                    file.WriteLine("{0},{1}", entry.Key, entry.Value);
                }
            }
        }

        /// <summary>
        /// Get all available options for the given graph.
        /// </summary>
        /// <param name="cand"></param>
        /// <returns></returns>
        public List<option> GetAvailableOptions(candidate cand) {
            var options = new List<option>();
            options.AddRange(Settings.rulesets[0].recognize(cand.graph, false));//ruleset[0] is non-carboxy ruleset
//            for (var i = options.Count - 1; i >= 0; i--) // trim options that would lead to invalid states
////                if (options[i].ruleNumber == CarboxylRule || !IsValidChild(cand, options[i])) 
//                if (options[i].ruleNumber == CarboxylRule) {
//                    options.RemoveAt(i);
//                    Console.WriteLine("Find carboxy rule in non-carbox ruleset!");
//                }
                //MC:IsValidChild actually creates child. So, it would seem you wouldn't
                // need to return options here, but rather all the valid children. If branching factor 
                    // is too large, then need to re-think this.
                 // Do not consider CarboxyOpt here
       
            var set = new HashSet<String>();
            for (var i = options.Count - 1; i >= 0; i--) {
                var evalCand = CopyAndApplyOption(options[i], cand, false);
                var smile = OBFunctions.moltoSMILES(OBFunctions.designgraphtomol(evalCand.graph));
                if (set.Contains(smile))
                {
                    options.RemoveAt(i);
                } else
                {
                    set.Add(smile);
                }
            }    
            return options;
        }
        
        /// <summary>
        /// Returns whether the graph violates angle or carboxyl-blocking constraints.
        /// </summary>
        private static bool IsValidChild(candidate cand, option opt) {
            var newCand = cand.copy();
            var newOpt = opt.copy();
            SearchProcess.transferLmappingToChild(newCand.graph, cand.graph, newOpt);
            newOpt.apply(newCand.graph, null);
            
            var mol = OBFunctions.designgraphtomol(newCand.graph);
            var mapping = OBFunctions.findcarboxylates(mol);
            if (mapping.Count < 2) return true; 
            
            var carbA = mol.GetAtom(mapping[0][1]); // carbon in carboxylate
            var aA = mol.GetAtom(mapping[0][3]); // atom that the carbon connects to
            var carbB = mol.GetAtom(mapping[1][1]);
            var aB = mol.GetAtom(mapping[1][3]);
            var pAngle = OBFunctions.pairwiseangle(carbA, aA, carbB, aB); // angle between carboxylates

            return !OBFunctions.carboxylatesBlocked(mol, aA, carbA, aB, carbB) && pAngle >= AngleFloor;

        }

        /// <summary>
        /// Returns the option corresponding to the carboxyl rule application.
        /// </summary>
        /// <param name="cand"></param>
        /// <returns></returns>
        public List<option> GetCarboxylOptions(candidate cand) {
//            cand.graph = Minimize(cand.graph);
            var options = new List<option>();
            options.AddRange(Settings.rulesets[1].recognize(cand.graph, false));//ruleset[1] is carboxy ruleset
//            for (var i = options.Count - 1; i >= 0; i--) // trim options that would lead to invalid states
//                if (options[i].ruleNumber != 1) {
//                    Console.WriteLine("Find non-carboxy rule {0} in carbox ruleset!", options[i].ruleNumber);
//                    options.RemoveAt(i);
//                }
                    // Do not consider CarboxyOpt here
            var set = new HashSet<String>();
            for (var i = options.Count - 1; i >= 0; i--) {
                var evalCand = CopyAndApplyOption(options[i], cand, false);
                var smile = OBFunctions.moltoSMILES(OBFunctions.designgraphtomol(evalCand.graph));
                if (set.Contains(smile))
                {
                    options.RemoveAt(i);
                } else
                {
                    set.Add(smile);
                }
            }
                    
            return options;// return the first one by default
        }

        /// <summary>
        /// Apply the option to the candidate and store the agent's evaluation.
        /// </summary>
        public void ApplyOption(option opt, candidate cand, bool doMinimize, bool doEvaluate=false) {
            cand.graph.globalVariables.Add(cand.f0); // track fitness values of previous states 
            
            // Minimize and evaluate
            opt.apply(cand.graph, null);
            cand.addToRecipe(opt);
            if (doMinimize)
                cand.graph = Minimize(cand.graph);
            if (doEvaluate)
                cand.f0 = Evaluate(cand);
            // MC Is this like line 64? by the way, tracking fitness here is better than in cand.graph
        }

        /// <summary>
        /// Copies the candidate graph, transfers the L-mapping, and returns the resultant candidate.
        /// </summary>
        public candidate CopyAndApplyOption(option opt, candidate cand, bool doMinimize) {
            var newCand = cand.copy();
            var newOpt = opt.copy();
            SearchProcess.transferLmappingToChild(newCand.graph, cand.graph, newOpt);
            ApplyOption(newOpt, newCand, doMinimize);
            return newCand;
        }
        
        /// <summary>
        /// Clean way to minimize a graph.
        /// </summary>
        private static designGraph Minimize(designGraph graph) {
            var mol = QuickMinimization(OBFunctions.designgraphtomol(graph), IODir + "rank" + ".lmpdat",
                IODir + "rank" + ".coeff", false, 0);
            OBFunctions.updatepositions(graph, mol);
            return OBFunctions.tagconvexhullpoints(graph);
        }

        private static OBMol QuickMinimization(OBMol mol, string coeff, string lmpdat, bool periodic, int rankMe) {
            double padding = 50;
            const double etol = 0.0;
            const double ftol = 1.0e-6;
            const int maxiter = 40000;
            const int maxeval = 20000;
            

            var minSettings = new lammps.LAMMPSsettings();
            if (periodic) {
                minSettings.boundary = "p p p";
                padding = 0;
            }

            if (File.Exists(coeff))
                File.Delete(coeff);
            if (File.Exists(lmpdat))
                File.Delete(lmpdat);
            Converter.moltoUFF(mol, coeff, lmpdat, false, padding);
 

            string[] lmparg = {"", "-screen", "none", "-log", "log.lammps." + rankMe};
            using (var lmps = new lammps(minSettings, lmparg)) {
                lmps.runCommand("read_data " + lmpdat);
                lmps.openFile(coeff);
                lmps.minimize(etol, ftol, maxiter, maxeval, "cg");
                OBFunctions.updatexyz(mol, lmps.getAtomPos());
            }


            return mol;
        }

        /// <summary>
        /// Evaluate the fitness of the molecule.
        /// </summary>
        /// <param name="cand"></param>
        /// <returns></returns>
        public double Evaluate(candidate cand) {
            return 0;
        }
        
        /// <summary>
        /// Check if a molecular has been evaluated or not.
        /// </summary>
        /// <param name="cand"></param>
        /// <returns></returns>
        public bool HasEvaluate(candidate cand) {
            return previous.ContainsKey(OBFunctions.moltoSMILES(OBFunctions.designgraphtomol(cand.graph)));
        }



        

        
        /// <summary>
        /// Calculate the angle between 2 carboxylates
        /// </summary>
        public static double CalAngle(OBMol mol) {
            
            
            var mapping = OBFunctions.findcarboxylates(mol);
            
            var carbA = mol.GetAtom(mapping[0][1]); // carbon in carboxylate
            var aA = mol.GetAtom(mapping[0][3]); // atom that the carbon connects to
            var carbB = mol.GetAtom(mapping[1][1]);
            var aB = mol.GetAtom(mapping[1][3]);
            var pAngle = OBFunctions.pairwiseangle(carbA, aA, carbB, aB); // angle between carboxylates

            return pAngle;
        }
        
        /// <summary>
        /// Whether the candidate has enough carboxylates and conforms to angle- and carboxylate-blocking constraints.
        /// </summary>
        public static bool CanCalculateReward(OBMol mol) {
            
            
            var mapping = OBFunctions.findcarboxylates(mol);
            if (mapping.Count < 2) return false; 
            
            var carbA = mol.GetAtom(mapping[0][1]); // carbon in carboxylate
            var aA = mol.GetAtom(mapping[0][3]); // atom that the carbon connects to
            var carbB = mol.GetAtom(mapping[1][1]);
            var aB = mol.GetAtom(mapping[1][3]);
            var pAngle = OBFunctions.pairwiseangle(carbA, aA, carbB, aB); // angle between carboxylates

            return !OBFunctions.carboxylatesBlocked(mol, aA, carbA, aB, carbB) && pAngle >= AngleFloor;
        }
        
        private void RunLammpsUsingSGE(string file, string extraArgs, string workingdir, int linknum) {
            var flagFile = Path.Combine(workingdir, "linker" + linknum + "-ave-force.d"); // says if lammps is finished
            if (File.Exists(flagFile)) {
                File.Delete(flagFile);
            }
            

            // Tell qsub to run with binary 
            const string lammpsExecutable = "/home/manion/expansion/Downloads/lammps-16Feb16/src/lmp_ubuntu";
            //var lammpsExecutable = Path.Combine(Settings.OutputDirAbs, "lmp_serial");
            var arguments = " -in " + file + " " + extraArgs;
            var fullArguments = lammpsExecutable + arguments;
            var jobname = "linker" + linknum + " -b y -wd " + workingdir + " "; //execute as binary
            OBFunctions.submitlammps(fullArguments, workingdir, jobname); // lmp -v Linker "linker-" -v LinkerRoot
                
            while (!File.Exists(flagFile)) {}  
//            Console.WriteLine(flagFile);
        }
    }
}