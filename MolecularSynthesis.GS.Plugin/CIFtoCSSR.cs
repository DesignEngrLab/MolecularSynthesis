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
using System.Timers;
using System.Threading.Tasks;
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
using System.Timers;
using System.Diagnostics;
using System.Threading.Tasks;
using System.Text;
using System;
using System.Diagnostics;
using System.ComponentModel;

namespace TestOpenBabel
{
    public class CIFtoCSSR : SearchProcess
    {
        public override string text => "CIFtoCSSR";
        static readonly Dictionary<string, int> elementTabel = new Dictionary<string, int>()
        {
            {"H", 1},
            {"C", 6},
            {"N", 7},
            {"O", 8},
            {"Cl", 17},
            {"Br", 35},
            {"F", 9}
        };
        //deault constructor
        public CIFtoCSSR(GlobalSettings settings) : base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 2;
            AutoPlay = true;
        }
        protected override void Run()
        {
            //var resultMol = OBFunctions.designgraphtomol(seedGraph);
            // after convert from designgraph to .mol file
            // save the result as .mol file
            // call minimize.exe to do the energy minimization
            // 
            var stopwatch = new Stopwatch();


            List<string> list = new List<string>();

            //string text = System.IO.File.ReadAllText(@"C:\Users\zhang\source\repos\tobacco_3.0\output_cifs\pcu_v1-6c_Zn_1_Ch_1-Candidate19_fer_tobacco.cif");
            //Console.WriteLine("-------------------------------------------------");


            // find the cell length
            string[] lines = System.IO.File.ReadAllLines(@"C:\Users\zhang\source\repos\tobacco_3.0\output_cifs\pcu_v1-6c_Zn_1_Ch_1-Candidate9_fer_tobacco.cif");
            //Console.WriteLine("----" +lines[9].Split(' ')[lines[9].Split(' ').Length-1] + "-----");
            //Console.WriteLine(lines[9].Split(' ')[20]);
            string[] words = lines[9].Split(' ');
            //-----------------------------------------------------------

            // writing the beginning part of the file
            list.Add("                  " + lines[9].Split(' ')[20] + "  " + lines[9].Split(' ')[20] + "  " + lines[9].Split(' ')[20]);
            list.Add("		90  90  90  SPGR =  1 P 1		 OPT = 1");


            //find total number of nodes, start from line 22
            int N = 22;

            while (lines[N].Contains("loop_") == false)
            {
                N = N + 1;
            }
            list.Add((N - 22).ToString() + "  " + "0");
            list.Add("0" + "      :");
            //  writing node information, 7 element, 10 X value, 13 Y value, 17 Z value
            N = 22;
            while (lines[N].Contains("loop_") == false)
            {
                List<string> cssr = new List<string>();
                foreach (var item in lines[N].Split(" "))
                {
                    if (item != string.Empty)
                    {
                        cssr.Add(item);
                        Console.WriteLine(item);
                    }
                }
                list.Add((N - 21).ToString() + " " + cssr[1] + " " + cssr[2] + " " + cssr[3] + " " + cssr[4]);
                //list.Add((N-21).ToString()+" "+ lines[N].Split(" ")[7]+ " "+ lines[N].Split(" ")[10] + " "+ lines[N].Split(" ")[13]+" "+ lines[N].Split(" ")[17]);

                N = N + 1;
            }
            System.IO.File.WriteAllLines(@"C:/Users/zhang/Desktop/testCifToCssr.cssr", list);



            //foreach (string line in list)
            //{
            //    // Use a tab to indent each line of the file.
            //    Console.WriteLine("\t" + line);
            //}

            //string result = string.Join(".", list);
            //Console.WriteLine("--------------------------------");

            //foreach (var item in nodeatomlookup)
            //{
            //    Console.WriteLine(item);
            //}
            //Console.WriteLine("--------------------------------");


            //Console.WriteLine($"RESULT: {result}");
            //File.WriteAllText("/nfs/hpc/share/zhangho2/MolecularSynthesis/examples/HS.txt", list);
            System.IO.File.WriteAllLines(@"/nfs/hpc/share/zhangho2/MolecularSynthesis/examples/HS__HS.cif", list);
            //------------------------------------------------------------------------------------



            //resultMol = OBFunctions.InterStepMinimize(resultMol);
            //OBFunctions.updatepositions(StartState.graph, resultMol);




            //var stopwatch = new Stopwatch();

            ////C:\Users\zhang\Desktop
            //string text = System.IO.File.ReadAllText(@"C:\Users\zhang\Desktop\PPP.cif");

            //Console.WriteLine("-------------------------------------------------");
            //System.Console.WriteLine("Contents of WriteText.txt = {0}", text);

            //string[] lines = System.IO.File.ReadAllLines(@"C:\Users\zhang\Desktop\PPP.cif");

            //// Display the file contents by using a foreach loop.
            //System.Console.WriteLine("Contents of WriteLines2.txt = ");
            //foreach (string line in lines)
            //{
            //    // Use a tab to indent each line of the file.
            //    Console.WriteLine("\t" + line);
            //}

















            //var resultMol = OBFunctions.designgraphtomol(seedGraph);
            //resultMol = OBFunctions.InterStepMinimize(resultMol);
            //OBFunctions.updatepositions(seedGraph, resultMol);

            //SearchIO.addAndShowGraphWindow(seedGraph);

            //var result = Evaluation.FindLengthAndRadius(seedGraph);
            ////[0] is Length, [1] is Radius, unit is 
            //SearchIO.output("Length is: " + result[0]);
            //SearchIO.output("Radius is: " + result[1]);


            //resultMol = OBFunctions.InterStepMinimize(resultMol);
            //OBFunctions.updatepositions(seedGraph, resultMol);
            //var FinalResultMol = OBFunctions.designgraphtomol(seedGraph);

            //var conv = new OBConversion();
            //conv.SetInAndOutFormats("pdb", "mol");
            //conv.WriteFile(FinalResultMol, Path.Combine("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\output", "Test102.mol"));



            //SearchIO.addAndShowGraphWindow(seedGraph);

            //var result = new double[2];
            //result=Evaluation.FindLengthAndRadius(seedGraph);
            ////[0] is Length, [1] is Radius, unit is 
            //SearchIO.output("Length is: "+ result[0]);
            //SearchIO.output("Radius is: "+ result[1]);

            //string[] PathOfXyz = Directory.GetFiles(@"C:\Users\zhang\desktop", "*.xyz");
            //string CifFilename = Path.GetFileName(PathOfXyz[0]);
            //foreach (var namexxx in PathOfXyz)
            //{
            //    Console.WriteLine(Path.GetFileName(namexxx));
            //}
            //Console.WriteLine(Path.GetFileName(PathOfXyz[0]));
            //Console.WriteLine(Path.GetFileName(PathOfXyz[1]));

            //RunBat("C:\\Users\\zhang\\source\\repos\\PoreBlazer\\Windows\\HKUST1\\run.bat");

            //C:\Users\zhang\source\repos\PoreBlazer\Windows\HKUST1\results.txt
            //string[] EvaluationResult = System.IO.File.ReadAllLines(@"C:\Users\zhang\source\repos\PoreBlazer\Windows\HKUST1\results.txt");
            //Console.WriteLine(EvaluationResult[174]);

            stopwatch.Restart();
            var elapsed = stopwatch.Elapsed;
            Console.WriteLine("completed in {0}", elapsed);

            //var conv = new OBConversion();
            //conv.SetInAndOutFormats("pdb", extension);
            //conv.WriteFile(resultMol, filename + "." + extension);
            //File.AppendText()
        }

        private void RunBat(string batPath)
        {
            Process pro = new Process();

            FileInfo file = new FileInfo(batPath);
            pro.StartInfo.WorkingDirectory = file.Directory.FullName;
            pro.StartInfo.FileName = batPath;
            pro.StartInfo.CreateNoWindow = false;
            pro.StartInfo.UseShellExecute = true;
            pro.Start();
            pro.WaitForExit();
        }
    }




}

