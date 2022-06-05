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
    public class random_search : SearchProcess
    {
        // give desiredMoment
        // [] desiredMoment = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        //static double[] desiredLenghtAndRadius = new double[] { 300, 50 };
        static Random rnd = new Random(11);
        //static double[] desiredLenghtAndRadius = new double[] { 245.277, 89.53 };
        // RS0 R7 R3; RS1 R1 R2

        static double[] desiredLenghtAndRadius = new double[] { 560, 140 };
        // RS0 R3 R4 R5 R6

        public random_search(GlobalSettings settings) : base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 2;
            AutoPlay = true;
        }

        static TreeCandidate noneparallel = new TreeCandidate(new candidate());

        public override string text
        {
            get { return "random_search"; }
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
            int TotalNumber = 30;
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
            try
            {

                Parallel.For(0, TotalNumber, count =>
                {

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
                        var totalFGnumber = rand.Next(0, PossibleFGnumber/9);


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

                    var score = Evaluation.distance(candidate, desiredLenghtAndRadius);

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

                    var FinalResultMol = OBFunctions.designgraphtomol(candidate.graph);

                    var conv = new OBConversion();
                    lock (noneparallel)
                        conv.SetInAndOutFormats("pdb", "mol");

                    string name = ".mol";
                    name = "HaHa" + Convert.ToString(System.Threading.Thread.CurrentThread.ManagedThreadId) + name;
                    lock (noneparallel)
                        conv.WriteFile(FinalResultMol, Path.Combine("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\examples", name));

                    //}
                });

            }


            catch (Exception exc)
            {
                Console.WriteLine("exc.message");
            }

                        

            var resultMol = OBFunctions.designgraphtomol(bestCandidate.graph);
            resultMol = justMinimize(resultMol);
            OBFunctions.updatepositions(bestCandidate.graph, resultMol);

            var score = Evaluation.distance(bestCandidate, desiredLenghtAndRadius);

            Console.WriteLine("---------------------------------------");
            Console.WriteLine(score);
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

            System.IO.File.WriteAllLines(@"C:\Users\zhang\source\repos\MolecularSynthesis\examples\RandomSearchRecord.txt", resultCollector);


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
                conv.WriteFile(mol, Path.Combine("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\examples", filename));
            string minimizeOutput;

            using (Process proc = new Process())
            {
                // C:\Program Files (x86)\OpenBabel-3.1.1

                proc.StartInfo.FileName = "C:\\Program Files (x86)\\OpenBabel-3.1.1\\obminimize.exe";
                proc.StartInfo.Arguments = "-ff GAFF " + filename;

                //proc.StartInfo.Arguments = "-n200 minimize.mol"; //can add arguments here like number of iterations,
                // or '-c' convergence criteria
                proc.StartInfo.ErrorDialog = false;
                proc.StartInfo.WorkingDirectory = "C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\examples";
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






    }

}






