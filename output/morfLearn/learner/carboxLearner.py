from . import util
import numpy as np
import os
import random
import torch
from torch.utils.data import Dataset, DataLoader
import pickle
import sys



class CarboxLearner(object):
	def __init__(self, data_dir, task, model='point'):
		self.data_dir = data_dir
		self.task = task
		self.model = model
		assert model in ['point']
		self.valueNet = None
		if model == 'point':
			self.valueNet = util.point.PointNet(self.task)
		self.data_set = {}

	def predict(self, linkerName):
		arr = np.load(os.path.join(self.data_dir, "features", self.model, linkerName + ".npy"))
		return arr.shape


	def learn_fix(self, batch_size=200, epoch=50):
		def create_dataset(task):
			if self.task == "Classification":
				pkl_file = "pointBFSCsf.pkl"
				label_file = os.path.join(self.data_dir, "property", "stiffC", "label.npy")
			elif self.task == "Regression":
				pkl_file = "pointBFSReg.pkl"
				label_file = os.path.join(self.data_dir, "property", "stiff.npy")
				
			if pkl_file in os.listdir(os.path.join(self.data_dir, "features")):
				self.data_set = pickle.load(open(os.path.join(self.data_dir, "features", pkl_file), "rb"))
				print("BFS Dataset Loaded")
			else:
				self.data_set = {}
				data_label = np.load(label_file)
				point_dir = os.path.join(self.data_dir, "features", "pointBFS")
				for npArray in os.listdir(point_dir):
					key = int(npArray[npArray.find("pointCloud") + len("pointCloud") : npArray.find(".npy")])
					x = np.load(os.path.join(point_dir, npArray))
					y = data_label[key]
					self.data_set[key] = (x,y)
				pickle.dump(self.data_set, open(os.path.join(self.data_dir, "features", pkl_file), "wb"))
				print("BFS Dataset Created")
				
		def save(task, net):
			if task == "Classification":
				model_file = os.path.join("model", "carbox_value_" + self.model + "_csf.pt")
			elif task == "Regression":
				model_file = os.path.join("model", "carbox_value_" + self.model + "_reg.pt")
			if "model" not in os.listdir(os.getcwd()):
				os.mkdir("model")
			torch.save(net.state_dict(), model_file)
			print("Save Model")
	
		def load(task, net):
			if task == "Classification":
				model_file = os.path.join("model", "carbox_value_" + self.model + "_csf.pt")
			elif task == "Regression":
				model_file = os.path.join("model", "carbox_value_" + self.model + "_reg.pt")
			net.load_state_dict(torch.load(model_file))
			print("Model Loaded")
			
		def evaluate(data_loader, task, net):
			sum_loss = 0
			sum_acc = 0
			n_batch = 0
			for i,(x,y) in enumerate(data_loader):
				y_hat = net(x)
				loss = net.criterion(y_hat, y).detach().cpu().numpy()
				sum_loss += loss
				if task == "Classification":
					acc = torch.mean((torch.argmax(y_hat, dim=1) == torch.argmax(y, dim=1)).float()).detach().cpu().numpy()
					sum_acc += acc
				n_batch = i+1
			if task == "Classification":
				return sum_loss/n_batch, sum_acc/n_batch
			return sum_loss/n_batch
				
			
		create_dataset(self.task)
		# load(self.task, self.valueNet)
		train_idx = np.load(os.path.join(self.data_dir, "trainIdxNoTrim.npy"))
		valid_idx = np.load(os.path.join(self.data_dir, "testIdxNoTrim.npy"))
		sample = list(map(lambda i: self.data_set[i], np.arange(len(self.data_set))))
		x = np.array(list(map(lambda x: x[0], sample)), dtype=np.float32)
		y = np.array(list(map(lambda x: x[1], sample)), dtype=np.float32)
		train_set = util.point.Point_dataset(x[train_idx], y[train_idx])
		valid_set = util.point.Point_dataset(x[valid_idx], y[valid_idx])
		train_loader = DataLoader(dataset=train_set, batch_size=batch_size, shuffle=True)
		eval_loader = DataLoader(dataset=valid_set, batch_size=2139, shuffle=False)
		
		if self.task == "Classification":
			best_validation = 0
		else:
			best_validation = float('Inf')
		for e in range(epoch):
			for i,(x,y) in enumerate(train_loader):
				y_hat = self.valueNet(x)
				loss = self.valueNet.criterion(y_hat, y)
				self.valueNet.optimizer.zero_grad()
				loss.backward()
				self.valueNet.optimizer.step()
				
				
			if self.task == "Classification":
				train_loss, train_acc = evaluate(train_loader, self.task)
				valid_loss, valid_acc = evaluate(eval_loader, self.task)
				print("Epoch %d: training loss %.3f, training acc %.3f, validation loss %.3f, validation acc %.3f." 
			% (e, train_loss, train_acc, valid_loss, valid_acc))
				if valid_acc > best_validation:
					save(self.task, self.valueNet)
					best_validation = valid_acc
			elif self.task == "Regression":
				train_loss = evaluate(train_loader, self.task, self.valueNet)
				valid_loss = evaluate(eval_loader, self.task, self.valueNet)
				print("Epoch %d: training loss %.3f, validation loss %.3f." 
			% (e, train_loss, valid_loss))
				if valid_loss < best_validation:
					save(self.task, self.valueNet)
					best_validation = valid_loss
			

			sys.stdout.flush()
			
		







