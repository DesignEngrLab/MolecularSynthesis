using OpenBabel;
using System;
using OpenBabelFunctions;
using GraphSynth.Search;
using GraphSynth;
using MolecularSynthesis.Plugin;
using System.Xml.XPath;


namespace TestOpenBabel
{
    public class TestBenzRing : SearchProcess
    {
        public override string text => "TestBenz";

        //deault constructor
        public TestBenzRing(GlobalSettings settings) : base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 2;
            AutoPlay = true;
        }
        protected override void Run()
        {
            var resultMol= OBFunctions.designgraphtomol(seedGraph);
            resultMol = OBFunctions.InterStepMinimize(resultMol);
            OBFunctions.updatepositions(seedGraph, resultMol);
            SearchIO.addAndShowGraphWindow(seedGraph);

            var result = new double[2];
            result=Evaluation.FindLengthAndRadius(seedGraph);
            //[0] is Length, [1] is Radius
            SearchIO.output(result[0]+ " " +result[1]);

            var conv = new OBConversion();
            conv.SetInAndOutFormats("pdb", "mol");
            conv.WriteFile(resultMol, @"C:\Users\zhang\source\repos\MolecularSynthesis\output\minimize.mol");
        }
    }




}

