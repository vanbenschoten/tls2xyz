#cctbx.python tls2xyz_run.py pdb=XXX.pdb models=10
from iotbx import pdb
import sys
import math
import mmtbx.tls.tools
import tls_as_xyz
from cctbx.array_family import flex
from cctbx import crystal
from cctbx import miller
from iotbx import scalepack
import os

def run(arguments):

	args = get_input_dict(arguments)

	file = input_model(args['pdb'])
	file.get_tls_ensembles(int(args['models']))
	file.combine_tls_models()



def get_input_dict(args):
	dic = dict()
	for arg in args:
		spl=arg.split('=')
		if len(spl)==2:
			dic[spl[0]] = spl[1]

	return dic


class input_model:

	def __init__(self, pdb):
		self.pdb = pdb
		self.extract_tls_information()

	def extract_tls_information(self):

		pdb_inp = pdb.input(file_name=self.pdb)
  		h = pdb_inp.construct_hierarchy()

  		self.unit_cell = pdb_inp.crystal_symmetry_from_cryst1()

  		pdb_obj = pdb.hierarchy.input(file_name=self.pdb)
  		remark_3_records = pdb_inp.extract_remark_iii_records(3)
  		tls_extract = mmtbx.tls.tools.tls_from_pdb_inp(
    		remark_3_records = remark_3_records, 
    		pdb_hierarchy = h)

  		self.tls_segments = list()
  		self.tls_objects = list()
  		for i, tls_object in enumerate(tls_extract.tls_params):
  			self.tls_segments.append(tls_object.selection_string)
  			self.tls_objects.append(tls_object)

  		count = 0

  		self.atom_files = list()

		for group in self.tls_segments:

			sel = h.atom_selection_cache().selection(string = group)
			h_new = h.select(sel)

			#Write atomic coordinates of each 
			h_new.write_pdb_file(file_name='model_%d.pdb' %count)

			self.atom_files.append('model_%d.pdb' %count)

			count += 1

	def get_tls_ensembles(self,models):

		length = len(self.atom_files)

		self.ensemble_files = list()

		for i in range(0,length):
			tls_as_xyz.run(pdb_file_name = self.atom_files[i], tls_object=self.tls_objects[i], crystal_info=self.unit_cell, n_models=models, output_file_name='model_%d_tls.pdb' %i)
			self.ensemble_files.append('model_%d_tls.pdb' %i)

	def combine_tls_models(self):
		data = dict()

		data['original'] = dict()
		data['original']['pdb object'] = pdb.input(file_name=self.pdb)
		data['original']['hierarchy'] = data['original']['pdb object'].construct_hierarchy()
		data['original']['models'] = []
		data['original']['chains'] = []

		for model in data['original']['hierarchy'].models():
			data['original']['models'].append(model.id)

		for chain in data['original']['hierarchy'].chains():
			data['original']['chains'].append(chain.id)


		chains = []
		for item in data['original']['chains']:
			if item not in chains:
				chains.append(item)

		#Populating the dictionary with necessary info about each TLS model
		for model in self.ensemble_files:

			data[model] = dict()
			data[model]['pdb object'] = pdb.input(file_name=model)
			data[model]['hierarchy'] = data[model]['pdb object'].construct_hierarchy()
			data[model]['models'] = []

			count = 0
			for model_h in data[model]['hierarchy'].models():
				data[model]['models'].append(model_h)
				count += 1

			data['Number of models'] = count

		#Define the new PDB hierarchy to be written to
		output = pdb.hierarchy.ext.root()

		m_count = 0

		for i in range(0,data['Number of models']):
			new_model = pdb.hierarchy.ext.model()
			new_model.id = str(m_count)

		#Populate the new model with the chains present in the original model
			chain_list = []
			for item in chains:
				new_chain = pdb.hierarchy.ext.chain()
				new_chain.id = item
				chain_list.append(new_chain)

			for model in self.ensemble_files:
				m_working = data[model]['models'][m_count]
				m_chains = m_working.chains()

				for ch in m_chains:
					name = ch.id

					for nm_c in chain_list:
						nm_name = nm_c.id

						if name == nm_name:
							resi = ch.residue_groups()

							for r in resi:
								re = r.detached_copy()
								nm_c.append_residue_group(re)

		#Once all the atoms across all the models are appended and accounted for, attach all the chains to the new model
			for nm_c in chain_list:
				new_model.append_chain(nm_c)

			output.append_model(new_model)
			m_count += 1



		output.write_pdb_file(file_name='output.pdb')

		fin = open('output.pdb', 'r')
		fout = open('combined_models.pdb', 'w')

		lines = fin.readlines()

		for line in lines:
			data = line.split()
			if data[0] != 'BREAK':
				fout.write(line)

		fin.close()
		fout.close()
		os.remove('output.pdb')



if __name__ == '__main__':
	import sys
	run(arguments=sys.argv[1:])