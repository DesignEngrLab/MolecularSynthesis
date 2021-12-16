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
using System.Diagnostics;
using System;
using System.Threading;
using System.Threading.Tasks;

namespace MolecularSynthesis.GS.Plugin
{
    public class MCTSSimpleCampbell : SearchProcess
    {
        // give desiredMoment
        // [] desiredMoment = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        // RS0 R7 R3; RS1 R1 R2
        //RS0 R3 L=134;
        //static double[] desiredLenghtAndRadius = new double[] { 245.277, 89.53 };
        static Random rnd = new Random();
        candidate bestCandidate = null;
        double ExploreWeight = 50;
        static double[] desiredLenghtAndRadius = new double[] { 565, 140 };
        static double desiredPD = 25.07296;

        // RS0 R3 R4 R5 R6
        static TreeCandidate noneparallel = new TreeCandidate(new candidate());

        public MCTSSimpleCampbell(GlobalSettings settings) : base(settings)
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
            get { return "MCTSSimpleCampbell"; }
        }
        protected override void Run()
        {
            //var candidates = new SimplePriorityQueue<candidate, double>();

            // generate a random number 0 or 1 to decide the next rule is from RS0 or RS1
            //Random rnd = new Random();
            //rnd.Next(0, 2); // generate 0 or 1

            // use 10000 is that DS use 3000-70000 iteration for 9*9 go play , so guess 10000 is enough
            int iteration = 3000;

            //TreeCandidate node1 = new TreeCandidate() { S = 0, n=0, UCB=0 };

            // 1. check if this is the leaf node, if not, go to step 2 until it is a leaf node,if yes go to step 3
            // 2. find the children who has the best UCB value
            // 3. do random simulation
            // 4. update S,n,UCB value for the whole tree


            // rollout process should match linker length , add a collector to collect solutions from trees with different depth
            // change for loop to while loop, set stop criteria, like n=50
            // careful for result from evaluation, should include posive value and negative value


            int IterationTimes = 0;
            List<string> MCTSProcess = new List<string>();
            //List<string> resultCollector = new List<string>();



            var timer = new Stopwatch();
            timer.Start();

            //TreeCandidate current = StartState;

            //while (current.n<50)

            //Parallel.For(0, 1, count =>
            //{
            //Console.WriteLine($"value of count = {count}, thread = {Thread.CurrentThread.ManagedThreadId}");
            //Sleep the loop for 10 miliseconds
            //Thread.Sleep(10);

            TreeCandidate StartState = new TreeCandidate(seedCandidate);

            StartState.S = 0;
            StartState.n = 0;
            StartState.UCB = double.MaxValue;
            StartState.Children = new List<TreeCandidate>();
            StartState.Parent = null;

            List<string> resultCollector = new List<string>();
            double[] Everystep = new double[2];
            double score = 0;

            for (int i = 0; i < iteration; i++)
            {

                Console.WriteLine("-----------------------------iterationtime=" + i.ToString());
                // if abs(current.S - target value)  < stop criteria 
                //  record this recipe

                // need to save S value and n value, delete the added graph, back to StartState                                                  
                TreeCandidate current = StartState;
                /// SELECTION
                current = SelectPromisingNode(current);// until at leaf node 
                ///EXPANSION
                CreateAllChildrenStates(current);
                /// Rollout aka. Simulation
                if (current.n == 0)
                    (current.S, score) = Simulation(current);

                resultCollector.Add(score.ToString());

                //--------------------------------------------------------------------------------------
                // Back Propogation
                //--------------------------------------------------------------------------------------
                BackPropogation(FindAllParents(current), current);

                //IterationTimes = DisplayData(IterationTimes, MCTSProcess, current);



            }

            //save the best!!

            //ReportFinalData(StartState, MCTSProcess);
            var filename = "MCTS_data_forUse" + Thread.CurrentThread.ManagedThreadId.ToString();
            filename = filename + ".txt";
            System.IO.File.WriteAllLines(@"/nfs/hpc/share/zhangho2/MolecularSynthesis/examples/" + filename, resultCollector);

            //});

            timer.Stop();
            TimeSpan ts = timer.Elapsed;

            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
            ts.Hours, ts.Minutes, ts.Seconds,
            ts.Milliseconds / 10);
            Console.WriteLine("RunTime " + elapsedTime);



        }

        private void ReportFinalData(TreeCandidate StartState, List<string> MCTSProcess)
        {


            //TreeCandidate seed = new TreeCandidate(seedCandidate);
            var FinalResult = FinalRecipe(StartState);
            string SolutionIsDownBelow = "Solution is down below: ";
            MCTSProcess.Add(SolutionIsDownBelow);

            foreach (var option in FinalResult.recipe)
            {
                string SolutionInformation = option.ruleSetIndex + " " + option.ruleNumber + "------------";
                MCTSProcess.Add(SolutionInformation);
            }

            System.IO.File.WriteAllLines(@"C:\Users\zhang\source\repos\MolecularSynthesis\output\MCTSProcessRecord.txt", MCTSProcess);
        }

        private int DisplayData(int IterationTimes, List<string> MCTSProcess, TreeCandidate current)
        {
            IterationTimes = IterationTimes + 1;
            SearchIO.output("IterationTimes = ", IterationTimes);

            string times = "Iteration times: " + IterationTimes.ToString();
            MCTSProcess.Add(times);
            string CurrentSValue = "S = " + current.S.ToString();
            MCTSProcess.Add(CurrentSValue);
            string CurrentnValue = "n = " + current.n.ToString();
            MCTSProcess.Add(CurrentnValue);

            if (IterationTimes > 1)
            {

                string CurrentUCBValue = "UCB = " + CalculateUcb(current).ToString();
                MCTSProcess.Add(CurrentUCBValue);
            }
            string ChildrenNumber = "Children number = " + current.Children.Count.ToString() + "*************";
            MCTSProcess.Add(ChildrenNumber);
            string NumberOfRecipe = "number of recipe: " + current.recipe.Count.ToString();
            MCTSProcess.Add(NumberOfRecipe);
            string CurrentNodeRecipe = "current node recipe:";
            MCTSProcess.Add(CurrentNodeRecipe);

            if (current.recipe.Count == 0)
            {
                string NoRecipe = "no recipe" + "------------";
                MCTSProcess.Add(NoRecipe);
            }
            else
            {
                foreach (var option in current.recipe)
                {
                    //SearchIO.output(current.recipe[j].ruleSetIndex + " " + current.recipe[j].optionNumber);
                    string OptionInformation = option.ruleSetIndex.ToString() + " " + option.ruleNumber.ToString() + "------------";
                    MCTSProcess.Add(OptionInformation);
                }
            }

            string seperateline = "-------------------------------------------------------------------------";
            MCTSProcess.Add(seperateline);
            return IterationTimes;
        }

        public double CalculateUcb(TreeCandidate child)
        {
            if (child.n == 0)
                return double.MaxValue;
            else
                return child.S / child.n + ExploreWeight * Math.Sqrt(Math.Log(child.Parent.n) / child.n);
        }

        //create the bestchild as an intermidiate variable
        public TreeCandidate SelectPromisingNode(TreeCandidate current)
        {
            TreeCandidate bestChild = current;
            while (current.Children.Count != 0)
            {
                double bestUcb = double.MinValue;
                foreach (TreeCandidate child in current.Children.OrderBy(item => rnd.Next())) //random shuffle of children
                {
                    double Ucb = CalculateUcb(child);
                    child.UCB = Ucb;
                    if (Ucb > bestUcb)
                    {
                        bestUcb = Ucb;
                        bestChild = child;
                    }
                }
                current = bestChild;
            }
            return bestChild;
        }

        public TreeCandidate FinalRecipe(TreeCandidate StartState)
        {
            TreeCandidate bestChild = null;
            while (StartState.Children.Count != 0)
            {
                double bestS = double.MinValue;
                foreach (TreeCandidate child in StartState.Children)
                {

                    if (child.S / child.n > bestS)
                    {
                        bestS = child.S / child.n;
                        bestChild = child;
                    }
                }
                StartState = bestChild;

            }

            return bestChild;
        }


        public void CreateAllChildrenStates(TreeCandidate current)
        {
            // need to add one avaiable option from current ,add options into recipe
            var options = new[] { rulesets[0].recognize(current.graph), rulesets[1].recognize(current.graph) }.SelectMany(x => x);
            foreach (var opt in options)
            {
                var child = (TreeCandidate)current.copy();
                child.Parent = current;
                child.Children = new List<TreeCandidate>();
                child.n = 0;
                child.S = 0;
                child.UCB = double.MaxValue;
                transferLmappingToChild(child.graph, current.graph, opt);
                opt.apply(child.graph, null);
                child.addToRecipe(opt);
                current.Children.Add(child);
            }
        }


        public void BackPropogation(IEnumerable<TreeCandidate> parentpath, TreeCandidate current)
        {
            current.n++;
            foreach (TreeCandidate treeCandidate in parentpath)
            {
                treeCandidate.n++;
                treeCandidate.S += current.S;
            }
        }

        public (double, double) Simulation(TreeCandidate child)  //aka Simulation
        {
            OBMol resultMol = OBFunctions.designgraphtomol(child.graph);
            resultMol = justMinimize(resultMol);
            OBFunctions.updatepositions(child.graph, resultMol);

            var score = loss(child, desiredPD);
            var S = 100000 - score;
            if (bestCandidate == null || ((TreeCandidate)bestCandidate).S < S)
                bestCandidate = child;
            return (S, score);
        }
        public IEnumerable<TreeCandidate> FindAllParents(TreeCandidate current)
        {
            while (current.Parent != null)
            {
                var parent = current.Parent;
                yield return parent;
                current = parent;
            }
        }

        public int GetHeight(TreeCandidate current)
        {
            int height = 0;


            TreeCandidate child = (TreeCandidate)current.copy();
            //child.Parent = current;
            //child.Parent = current.Parent;

            while (child.Parent != null)
            {
                height = height + 1;
                child = child.Parent;
            }

            return height;
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
                conv.WriteFile(mol, Path.Combine("/nfs/hpc/share/zhangho2/MolecularSynthesis/examples", filename));
            string minimizeOutput;

            using (Process proc = new Process())
            {
                // C:\Program Files (x86)\OpenBabel-3.1.1

                proc.StartInfo.FileName = "/usr/local/apps/openbabel/3.1.1/bin/obminimize";
                proc.StartInfo.Arguments = "-ff UFF " + filename;

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

                proc.StartInfo.Arguments = " -res " + "output_framework.cssr";
                //C: \Users\zhang\source\repos\MolecularSynthesis\output
                proc.StartInfo.WorkingDirectory = "/nfs/hpc/share/zhangho2/zeo++-0.3";
                //C:\\Users\\zhang\\Desktop
                proc.StartInfo.RedirectStandardOutput = true;
                proc.StartInfo.UseShellExecute = false;
                proc.Start();

                proc.WaitForExit();
            }

            // 7. read data from relative file


            string contents = File.ReadAllText("/nfs/hpc/share/zhangho2/zeo++-0.3/output_framework.res");
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


            return Math.Pow(Convert.ToDouble(words[4]) - desiredPD, 2);
        }

    }
}






