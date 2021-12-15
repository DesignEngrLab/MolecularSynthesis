using System.Collections.Generic;
using System.Threading;
using GraphSynth;
using GraphSynth.Representation;
using GraphSynth.Search;
using Priority_Queue;
using System;
using System.Collections;
using System.Security.Cryptography.X509Certificates;
using System.IO;
using MolecularSynthesis.GS.Plugin;
using System.Linq;
using OpenBabel;
using OpenBabelFunctions;
using MolecularSynthesis.GS.Plugin;
using System.Diagnostics;
using System.Threading.Tasks;
using System.Threading;

namespace MolecularSynthesis.GS.Plugin
{
    public class random_search_AV : SearchProcess
    {
        // give desiredMoment
        // [] desiredMoment = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        //static double[] desiredLenghtAndRadius = new double[] { 300, 50 };
        static Random rnd = new Random(11);
        //static double[] desiredLenghtAndRadius = new double[] { 245.277, 89.53 };
        // RS0 R7 R3; RS1 R1 R2

        static double[] desiredLenghtAndRadius = new double[] { 560, 140 };
        // RS0 R3 R4 R5 R6
        static double desiredAV = 3.98;

        public random_search_AV(GlobalSettings settings) : base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 2;
            AutoPlay = true;
        }

        static TreeCandidate noneparallel = new TreeCandidate(new candidate());

        public override string text
        {
            get { return "random_search_AV"; }
        }

        protected override void Run()
        {

            //TreeCandidate StartState = new TreeCandidate(seedCandidate);

            //StartState.S = 0;
            //StartState.n = 0;
            //StartState.UCB = double.MaxValue;
            //StartState.Children = new List<TreeCandidate>();

            //var option0 = rulesets[0].recognize(StartState.graph);
            //var option1 = rulesets[1].recognize(StartState.graph);
            //var option2 = rulesets[2].recognize(StartState.graph);

            ////option0[6].apply(StartState.graph, null);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[6].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[6]);

            // add a stopwatch to record time
            var timer = new Stopwatch();
            timer.Start();

            // Randomly generate .mol and .xyz files
            int TotalNumber = 3000;
            var rand = new Random();

            TreeCandidate StartState = new TreeCandidate(seedCandidate);

            TreeCandidate bestCandidate = new TreeCandidate(seedCandidate);
            bestCandidate.S = 10000000;

            List<string> resultCollector = new List<string>();

            // To randomly generate candidate
            // 1. how many backbone rule within[0,5], then randomly generate that number of backbone rule 
            // 2. how many functional group rule, then randomly generate that number of functional group rule
            // 3. calculate distance bewteen target as evaluate function
            // 4. get the distribution of the tree, get estimated average and estimated varience for MCTS 
            //for (int i = 0; i < TotalNumber; i++)
            //{
            //try
            //{
                for (int i = 0; i < TotalNumber; i++)
                {


                    //Parallel.For(0, TotalNumber, count =>
                    //{

                    //for (int i = 0; i < TotalNumber; i++)

                    var candidate = (TreeCandidate)StartState.copy();

                    var option0 = rulesets[0].recognize(candidate.graph);
                    var option1 = rulesets[1].recognize(candidate.graph);
                    var option2 = rulesets[2].recognize(candidate.graph);


                    // generate backbone rules
                    var BBnumber = rand.Next(1, 5);

                    for (int j = 0; j < BBnumber; j++)
                    {
                        option0 = rulesets[0].recognize(candidate.graph);
                        var totalOptionNumber = option0.Count();
                        var OptionNumber = rand.Next(0, totalOptionNumber);
                        option0[OptionNumber].apply(candidate.graph, null);
                        candidate.addToRecipe(option0[OptionNumber]);
                    }


                    // generate functional group rules, remember to check if option1 is zero

                    option1 = rulesets[1].recognize(candidate.graph);

                    if (option1.Count == 0)
                    {
                        option2 = rulesets[2].recognize(candidate.graph);
                        option2[0].apply(candidate.graph, null);
                    }

                    else
                    {
                        var PossibleFGnumber = option1.Count;

                        // check how many position are available
                        var totalFGnumber = rand.Next(0, PossibleFGnumber / 9);


                        for (int j = 0; j < totalFGnumber; j++)
                        {
                            var OptionNumber2 = rand.Next(0, PossibleFGnumber);


                            option1[OptionNumber2].apply(candidate.graph, null);
                            candidate.addToRecipe(option1[OptionNumber2]);

                            option1 = rulesets[1].recognize(candidate.graph);
                            PossibleFGnumber = option1.Count;
                        }

                        option2 = rulesets[2].recognize(candidate.graph);
                        option2[0].apply(candidate.graph, null);
                        //candidate.addToRecipe(option2[0]);

                    }


                    var resultMol = OBFunctions.designgraphtomol(candidate.graph);
                    resultMol = justMinimize(resultMol);
                    OBFunctions.updatepositions(candidate.graph, resultMol);

                    var score = loss(candidate, desiredAV);

                    if (score < bestCandidate.S)
                    {
                        bestCandidate = candidate;
                        bestCandidate.S = score;

                    }
                    resultCollector.Add(score.ToString());

                    //----------------------------------------
                    //var score = Evaluation.TotalAtomMass(candidate);
                    //resultCollector.Add(score.ToString());


                    // check the candidate in .mol

                    //var FinalResultMol = OBFunctions.designgraphtomol(candidate.graph);

                    //var conv = new OBConversion();
                    //lock (noneparallel)
                    //    conv.SetInAndOutFormats("pdb", "mol");

                    //string name = ".mol";
                    //name = "HaHa" + Convert.ToString(System.Threading.Thread.CurrentThread.ManagedThreadId) + name;
                    //lock (noneparallel)
                    //    conv.WriteFile(FinalResultMol, Path.Combine("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\examples", name));

                    //}
                    //});
                }
            //}


            //catch (Exception exc)
            //{
              //  Console.WriteLine("exc.message");
           // }



            var bestresultMol = OBFunctions.designgraphtomol(bestCandidate.graph);
            bestresultMol = justMinimize(bestresultMol);
            OBFunctions.updatepositions(bestCandidate.graph, bestresultMol);

            var bestscore = Evaluation.distance(bestCandidate, desiredLenghtAndRadius);

            Console.WriteLine("---------------------------------------");
            Console.WriteLine(bestscore);
            Console.WriteLine("---------------------------------------");

            Console.WriteLine("---------------------------------------");
            foreach (var option in bestCandidate.recipe)
            {
                string SolutionInformation = option.ruleSetIndex + " " + option.ruleNumber + "------------" + option.optionNumber;
                Console.WriteLine(SolutionInformation);
            }
            Console.WriteLine("*******" + bestCandidate.S + "********");
            Console.WriteLine("---------------------------------------");

            timer.Stop();
            TimeSpan ts = timer.Elapsed;
            //string foo = "Time taken: " + timeTaken.ToString(@"m\:ss\.fff");

            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
            ts.Hours, ts.Minutes, ts.Seconds,
            ts.Milliseconds / 10);
            Console.WriteLine("RunTime " + elapsedTime);

            System.IO.File.WriteAllLines(@"/nfs/hpc/share/zhangho2/MolecularSynthesis/examples/RS_AV.txt", resultCollector);


            //}
        }

        private static OBMol justMinimize(OBMol mol)
        {
            var conv = new OBConversion();
            lock (noneparallel)
                conv.SetInAndOutFormats("pdb", "mol");

            //lock (noneparallel) 

            int ThreadNumber = System.Threading.Thread.CurrentThread.ManagedThreadId;
            string filename = "Test" + ThreadNumber.ToString() + ".mol";

            lock (noneparallel)
                conv.WriteFile(mol, Path.Combine("/nfs/hpc/share/zhangho2/MolecularSynthesis/examples/", filename));
            string minimizeOutput;

            using (Process proc = new Process())
            {
                // C:\Program Files (x86)\OpenBabel-3.1.1

                proc.StartInfo.FileName = "/usr/local/apps/openbabel/3.1.1/bin/obminimize";
                proc.StartInfo.Arguments = "-ff GAFF " + filename;

                //proc.StartInfo.Arguments = "-n200 minimize.mol"; //can add arguments here like number of iterations,
                // or '-c' convergence criteria
                proc.StartInfo.ErrorDialog = false;
                proc.StartInfo.WorkingDirectory = "/nfs/hpc/share/zhangho2/MolecularSynthesis/examples";
                //proc.StartInfo.RedirectStandardError = true;
                //proc.StartInfo.UseShellExecute = false;
                proc.StartInfo.RedirectStandardOutput = true;
                //proc.StartInfo.RedirectStandardInput = false;

                Console.Write("starting OBMinimize...");
                proc.Start();

                minimizeOutput = proc.StandardOutput.ReadToEnd();
                proc.WaitForExit();

            }
            lock (noneparallel)
                conv.ReadString(mol, minimizeOutput);

            //File.Delete("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\examples\\"+ filename);




            return mol;

        }

        private static double loss(candidate child, double desiredPD)
        {
            var FinalResultMol = OBFunctions.designgraphtomol(child.graph);

            var conv = new OBConversion();
            conv.SetInAndOutFormats("pdb", "mol");

            // 1. generate .mol file, move the zeo++ folder
            //int i = 911;
            string uniqueName = Guid.NewGuid().ToString("D");
            string name = ".mol";
            name = uniqueName + name;
            conv.WriteFile(FinalResultMol, Path.Combine("/nfs/hpc/share/zhangho2/MolecularSynthesis/output", name));

            string position1 = "/nfs/hpc/share/zhangho2/MolecularSynthesis/output/" + name;
            string position2 = "/nfs/hpc/share/zhangho2/zeo++-0.3/" + name;
            System.IO.File.Move(position1, position2, true);
            Console.WriteLine("\n");
            //2. .mol to.xyz

            string name2 = ".xyz";
            name2 = uniqueName + name2;
            string position3 = "/nfs/hpc/share/zhangho2/zeo++-0.3/" + name2;

            using (Process proc = new Process())
            {
                //"C:\Program Files\OpenBabel-3.1.1\obabel.exe"
                proc.StartInfo.FileName = "/usr/local/apps/openbabel/3.1.1/bin/obabel";
                proc.StartInfo.Arguments = position2 + " -O " + position3;
                proc.StartInfo.WorkingDirectory = "/nfs/hpc/share/zhangho2/zeo++-0.3";
                //C:\Users\zhang\Desktop
                //proc.StartInfo.RedirectStandardError = true;
                //proc.StartInfo.UseShellExecute = false;
                proc.StartInfo.RedirectStandardOutput = true;
                //proc.StartInfo.RedirectStandardInput = false;

                Console.Write("starting Converting...");
                proc.Start();

                //minimizeOutput = proc.StandardOutput.ReadToEnd();
                proc.WaitForExit();
            }

            Console.WriteLine("\n");

            // 3. get rid of two carboxylate
            string name3 = uniqueName + "_XXX" + ".xyz";
            string position4 = "/nfs/hpc/share/zhangho2/zeo++-0.3/" + name3;

            using (Process proc = new Process())
            {
                proc.StartInfo.FileName = "/nfs/hpc/share/zhangho2/zeo++-0.3/molecule_to_abstract";
                proc.StartInfo.Arguments = position3 + " 0 " + position4;
                proc.StartInfo.WorkingDirectory = "/nfs/hpc/share/zhangho2/zeo++-0.3";
                proc.StartInfo.RedirectStandardOutput = true;
                proc.Start();
                Console.Write("starting removing...");
                proc.WaitForExit();

            }
            Console.WriteLine("\n");
            // 4. read the new xyz file and copy

            string[] lines = System.IO.File.ReadAllLines(@position4);
            System.IO.File.WriteAllLines(@"/nfs/hpc/share/zhangho2/zeo++-0.3/ForEvaluation.xyz", lines);
            Console.WriteLine("Finish writing new .xyz file");
            Console.WriteLine("\n");

            // 5.  build MOF


            // ./framework_builder nets/pcu.cgd 1 output 6c_Zn_1_Ch.xyz ForEvaluation.xyz
            using (Process proc = new Process())
            {
                proc.StartInfo.FileName = "/nfs/hpc/share/zhangho2/zeo++-0.3/framework_builder";
                proc.StartInfo.Arguments = "nets/pcu.cgd 1 output 6c_Zn_1_Ch.xyz ForEvaluation.xyz";
                proc.StartInfo.WorkingDirectory = "/nfs/hpc/share/zhangho2/zeo++-0.3";
                proc.StartInfo.RedirectStandardOutput = true;
                proc.StartInfo.UseShellExecute = false;
                proc.StartInfo.RedirectStandardError = true;
                proc.Start();
                Console.Write("starting building...");
                proc.WaitForExit();

            }
            //  5.1 need to change the output file name for multithread
            //File.Delete("/nfs/hpc/share/zhangho2/zeo++-0.3/IOP.cssr");
            //System.IO.File.Move("/nfs/hpc/share/zhangho2/zeo++-0.3/output_framework.cssr", "/nfs/hpc/share/zhangho2/zeo++-0.3/IOP.cssr");




            // 6. find accessible volume 
            using (Process proc = new Process())
            {

                //C: \Users\zhang\source\repos\zeo++-0.3\network
                // /nfs/hpc/share/zhangho2/MolecularSynthesis/zeo++-0.3
                proc.StartInfo.FileName = "/nfs/hpc/share/zhangho2/zeo++-0.3/network";
                //proc.StartInfo.Arguments = name + " -O " + name2;

                proc.StartInfo.Arguments = " -vol 1.2 1.2 50000 " + "output_framework.cssr";
                //C: \Users\zhang\source\repos\MolecularSynthesis\output
                proc.StartInfo.WorkingDirectory = "/nfs/hpc/share/zhangho2/zeo++-0.3";
                //C:\\Users\\zhang\\Desktop
                proc.StartInfo.RedirectStandardOutput = true;
                proc.StartInfo.UseShellExecute = false;
                proc.Start();

                proc.WaitForExit();
            }

            // 7. read data from relative file


            string contents = File.ReadAllText("/nfs/hpc/share/zhangho2/zeo++-0.3/output_framework.vol");
            string[] words = contents.Split(' ');
            Console.WriteLine(contents);

            Console.WriteLine("------------------");
            foreach (var word in words)
            {
                Console.WriteLine(word);
            }

            //Console.WriteLine("Accessible volume:" + words[15]);
            //Console.WriteLine("Accessible Volume Fraction:  " + words[13]);
            //Console.WriteLine("Poresize: ", words[5]);
            //PoreSizeValue = Convert.ToDouble(words[5]);
            //Console.WriteLine("Poresizevalue: ", words[5]);

            //Results.Add("Accessible volume: " + words[15] + "---" + ThreadNumber.ToString());
            //Results.Add("Accessible Volume Fraction: " + words[13] + "---" + ThreadNumber.ToString());

            File.Delete("/nfs/hpc/share/zhangho2/zeo++-0.3/IOP.cssr");


            return Math.Pow(Convert.ToDouble(words[15]) - desiredAV, 2);
        }

    }
}






