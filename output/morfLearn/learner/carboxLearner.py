from . import util
import numpy as np
import os
import random
import torch
from torch.utils.data import Dataset, DataLoader
import pickle
import sys

if torch.cuda.is_available():
	print("Using GPU")
	Tensor = torch.cuda.FloatTensor
else:
	print("Using CPU")
	Tensor = torch.FloatTensor

class CarboxLearner(object): 
	def __init__(self, data_dir, task, feature="point", property="stiff"):
		self.data_dir = data_dir
		self.task = task
		self.feature = feature
		self.property = property
		assert self.feature in ["point"]
		assert self.property in ["stiff"]
		self.valueNet = None
		if feature == "point":
			self.valueNet = util.point.PointNet(self.task)
		self.data_set = {}

	def predict(self, linkerName):
		feature_file = os.path.join(self.data_dir, "feature", self.feature, linkerName + ".npy")
		arr = np.load(feature_file)
		arr = np.expand_dims(arr, axis=0)
		est = self.valueNet(Tensor(arr))
		return est.detach().cpu().numpy().squeeze()

	def addData(self, linkerName):
		property_file = os.path.join(self.data_dir, "property", self.property, linkerName + ".npy")
		self.data_set[linkerName] = np.load(property_file)
		return self.data_set




	