using System;
using System.IO;
using GraphSynth.Representation;
using OpenBabelFunctions;

namespace GraphSynth.Search.Algorithms {

    /// <summary>
    /// Randomly select a valid option at any step.
    /// </summary>
    public class Random : AbstractAlgorithm
    {

        public Random(GlobalSettings settings) : base(settings)
        {
        }


        protected override string RunDirectoryName => "Random";

        public option ChooseOption(candidate cand)
        {
            var options = GetAvailableOptions(cand);
            return options.Count > 0 ? options[Rand.Next(options.Count)] : null;

        }

        public option ChooseCarboxOption(candidate cand)
        {
            var options = GetCarboxylOptions(cand);
            return options.Count > 0 ? options[Rand.Next(options.Count)] : null;

        }

        /// <summary>
        /// Build a candidate molecule by taking size actions.
        /// </summary>
        public void buildTrail(candidate seed, int size, int trails)
        {
            for (var t = 0; t < trails; t++)
            {
                Console.WriteLine("Trail: {0}", t);
                var cand = seed.copy();
                var linkerDir = _runDirectory + "/trail" + t + "/";
                Directory.CreateDirectory(linkerDir);
                var coeff = Path.Combine(linkerDir, "linker0.coeff");
                var lmpdat = Path.Combine(linkerDir, "linker0.lmpdat");
                // Set up UFF and run lammps
                Converter.moltoUFF(OBFunctions.designgraphtomol(cand.graph), coeff, lmpdat, false, 100);

                var candSMILE = OBFunctions.moltoSMILES(OBFunctions.designgraphtomol(cand.graph));
                
                Console.WriteLine("{0}", candSMILE);
                for (var step = 1; step < size + 1; step++)
                {
                    var opt = ChooseOption(cand);
                    if (opt == null)
                        return;
                    ApplyOption(opt, cand, true);
                    var mol = OBFunctions.designgraphtomol(cand.graph);
                    candSMILE = OBFunctions.moltoSMILES(mol);
                    Console.WriteLine("{0}", candSMILE);
                    coeff = Path.Combine(linkerDir, "linker" + step + ".coeff");
                    lmpdat = Path.Combine(linkerDir, "linker" + step + ".lmpdat");

                    Converter.moltoUFF(mol, coeff, lmpdat, false, 100);
                }
                var carboxOpt = ChooseCarboxOption(cand);
                if (carboxOpt == null)
                    return;
                ApplyOption(carboxOpt, cand, true);
                candSMILE = OBFunctions.moltoSMILES(OBFunctions.designgraphtomol(cand.graph));
                Console.WriteLine("{0}", candSMILE);
                coeff = Path.Combine(linkerDir, "linker" + (size+1) + ".coeff");
                lmpdat = Path.Combine(linkerDir, "linker" + (size+1) + ".lmpdat");

                Converter.moltoUFF(OBFunctions.designgraphtomol(cand.graph), coeff, lmpdat, false, 100);
            }

        }
    }



}