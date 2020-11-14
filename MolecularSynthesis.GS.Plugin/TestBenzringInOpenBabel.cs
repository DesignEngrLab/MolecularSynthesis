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

namespace TestOpenBabel
{
    public class TestBenzRing : SearchProcess
    {
        public override string text => "TestBenz";
        public string filename = @"..\..\..\..\ForCiftest";
        public string extension = "xyz";
        //deault constructor
        public TestBenzRing(GlobalSettings settings) : base(settings)
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

            TreeCandidate StartState = new TreeCandidate(seedCandidate);

            StartState.S = 0;
            StartState.n = 0;
            StartState.UCB = double.MaxValue;
            StartState.Children = new List<TreeCandidate>();

            var option0 = rulesets[0].recognize(StartState.graph);
            var option1 = rulesets[1].recognize(StartState.graph);
            var option2 = rulesets[2].recognize(StartState.graph);

            option0 = rulesets[0].recognize(StartState.graph);
            option0[0].apply(StartState.graph, null);
            StartState.addToRecipe(option0[0]);

            option0 = rulesets[0].recognize(StartState.graph);
            option0[5].apply(StartState.graph, null);
            StartState.addToRecipe(option0[5]);

            option0 = rulesets[0].recognize(StartState.graph);
            option0[5].apply(StartState.graph, null);
            StartState.addToRecipe(option0[5]);

            option0 = rulesets[0].recognize(StartState.graph);
            option0[5].apply(StartState.graph, null);
            StartState.addToRecipe(option0[5]);

            option0 = rulesets[0].recognize(StartState.graph);
            option0[5].apply(StartState.graph, null);
            StartState.addToRecipe(option0[5]);

            option2 = rulesets[2].recognize(StartState.graph);
            option2[0].apply(StartState.graph, null);
            StartState.addToRecipe(option2[0]);

            //option1 = rulesets[1].recognize(StartState.graph);
            //option1[8].apply(StartState.graph, null);
            //StartState.addToRecipe(option1[8]);

            //option1 = rulesets[1].recognize(StartState.graph);
            //option1[26].apply(StartState.graph, null);
            //StartState.addToRecipe(option1[26]);

            var resultMol = OBFunctions.designgraphtomol(StartState.graph);
            resultMol = OBFunctions.InterStepMinimize(resultMol);
            OBFunctions.updatepositions(StartState.graph, resultMol);

            List<string> list = new List<string>();
            list.Add("Hi");
            list.Add("XXX");
            System.IO.File.WriteAllLines(@"C:\Users\zhang\source\repos\MolecularSynthesis\output\MCTSProcess.txt", list);

            //resultMol = OBFunctions.InterStepMinimize(resultMol);
            //OBFunctions.updatepositions(seedGraph, resultMol);
            //var FinalResultMol = OBFunctions.designgraphtomol(seedGraph);

            //var conv = new OBConversion();
            //conv.SetInAndOutFormats("pdb", "mol");
            //conv.WriteFile(FinalResultMol, Path.Combine("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\output", "Test102.mol"));



            SearchIO.addAndShowGraphWindow(seedGraph);

            var result = new double[2];
            result=Evaluation.FindLengthAndRadius(seedGraph);
            //[0] is Length, [1] is Radius, unit is 
            SearchIO.output("Length is: "+ result[0]);
            SearchIO.output("Radius is: "+ result[1]);

            //var conv = new OBConversion();
            //conv.SetInAndOutFormats("pdb", extension);
            //conv.WriteFile(resultMol, filename + "." + extension);
            //File.AppendText()
        }
    }




}

