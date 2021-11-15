import os
import sys
import csv

import numpy as np
from MDAnalysis.analysis import align

import frame
import initialize as init

logger = init.get_logger(__file__)
np.seterr(all='ignore')

class FrameProcessor:
	_FLOAT_PRECISION = np.float16
	_DECIMALS = 3
	_SEPARATOR = '\t'
	_HEADERS = None

	def __init__(self, structure, selected_atoms, *args, **kwargs):
		self.structure = structure
		self.selected_atoms = selected_atoms

		self.float_precision = kwargs.get('float_precision', self._FLOAT_PRECISION)
		self.decimals = kwargs.get('decimals', self._DECIMALS)
		self.separator = kwargs.get('separator', self._SEPARATOR)
		self.group_name = kwargs.get('group_name', None)

		self.filepath = kwargs.get('filepath', None)

	def initialize(self):
		if self.filepath is not None:
			self.file = open(self.filepath, 'w', newline='\n', encoding='utf-8')
		else:
			self.file = open(os.devnull, 'w', newline='\n', encoding='utf-8')
		self.table = csv.writer(self.file, delimiter=self.separator)

		if self._HEADERS is not None:
			self.table.writerow(self._HEADERS)

	def do(self):
		raise NotImplementedError

	def write(self):
		raise NotImplementedError

	def finish(self):
		self.file.close()


class Alignment(FrameProcessor):
	_HEADERS = ('Frame', r'RMSD ($\AA$)')

	def __init__(self, structure, mobile_atoms, target_atoms, *args, **kwargs):
		super().__init__(structure, mobile_atoms, *args, **kwargs)
		self.target_atoms = target_atoms
		self.rmsds = np.zeros(self.structure.n_frames, dtype=self.float_precision)

	def do(self):
		try:
			R, rmsd = align.rotation_matrix(
				self.selected_atoms.positions - self.selected_atoms.center_of_geometry(),
				self.target_atoms.positions - self.target_atoms.center_of_geometry()
			)
		except ValueError:
			# logger.error(f'Error while superposing with {self.structure.label}')
			return None
		self.structure.translate(-self.selected_atoms.center_of_geometry())
		self.structure.rotate(R)
		self.structure.translate(self.target_atoms.center_of_geometry())        

		self.rmsds[self.structure.frame] = rmsd

	def write(self):
		self.table.writerow([
			self.structure.frame+1,
			np.around(self.rmsds[self.structure.frame], self.decimals)
		])

	def finish(self):
		self.table.writerow(['Average', np.around(self.rmsds.mean(), self.decimals)])
		self.table.writerow(['Standard deviation', np.around(self.rmsds.std(), self.decimals)])
		self.table.writerow(['Median', np.around(np.median(self.rmsds), self.decimals)])

		super().finish()


class FrameWriter(FrameProcessor):
	_TOPOLOGY_FORMATS = {
		frame.Static: 'pdb',
		frame.Multiple: 'pdb'
	}

	_COORDINATE_FORMATS = {
		frame.Static: 'pdb',
		frame.Multiple: 'ncdf'    
	}

	def __init__(self, structure, selected_atoms, traj_filepath, *args, **kwargs):
		return None
		super().__init__(structure, selected_atoms, *args, **kwargs)    
		self.traj_filepath = os.path.splitext(traj_filepath)[0]
		self.topology_format = self._TOPOLOGY_FORMATS[type(self.structure)]
		self.coordinate_format = self._COORDINATE_FORMATS[type(self.structure)]

		self.topology_filepath = f'{self.traj_filepath}.{self.topology_format}'
		self.coordinate_filepath = f'{self.traj_filepath}.{self.coordinate_format}'

		self.multiframe = isinstance(self.structure, frame.Multiple)

	def initialize(self):
		return None
		super().initialize()
		self.traj_file = MDAnalysis.Writer(self.coordinate_filepath, self.structure.protein.n_atoms, multiframe=self.multiframe)

		if self.multiframe:
			self.selected_atoms.write(self.topology_filepath)

	def do(self):
		return None
		self.traj_file.write(self.selected_atoms)

	def write(self):
		pass

	def finish(self):
		return None
		self.traj_file.close()
		super().finish()


class FramePropertyExtractor(FrameProcessor):
	def __init__(self, structure, selected_atoms, *args, **kwargs):
		super().__init__(structure, selected_atoms, *args, **kwargs)
		self.values = np.zeros(self.shape, dtype=self.float_precision)

	@property
	def shape(self):
		return self.structure.n_frames, self.select_atoms.n_atoms, 3
	

class FrameCoordinate(FramePropertyExtractor):
	def __init__(self, structure, selected_atoms, *args, **kwargs):
		super().__init__(structure, selected_atoms, *args, **kwargs)

		# Writing atomic coordinates is not implemented for the moment
		self.file = open(os.devnull, 'w')

	def do(self):
		self.values[self.structure.frame] = self.atoms.positions

	def write(self):
		pass


class FrameCentroid(FramePropertyExtractor):
	_HEADERS = ('Frame', 'X', 'Y', 'Z')

	def __init__(self, structure, selected_atoms, *args, **kwargs):
		super().__init__(structure, selected_atoms, *args, **kwargs)

	@property
	def shape(self):
		return self.structure.n_frames, 3

	def do(self):
		self.values[self.structure.frame] = self.selected_atoms.center_of_geometry()

	def write(self):
		self.table.writerow(
			[self.structure.frame+1] + 
			list(np.around(self.values[self.structure.frame], self.decimals)))

	def finish(self):
		self.table.writerow(['Average'] + list(np.around(self.values.mean(axis=0), self.decimals)))
		super().finish()


class ResultEntry:
	def __init__(self):
		self.structure = None
		self.values = np.array([])


class ProcessorResult:
	def __init__(self, *frame_processors):
		self.structures = []
		self.group_names = []
		self.values = np.array([])

		self.add(*frame_processors)

	@property
	def structure(self):
		if len(self.structures) == 1:
			return self.structures[0]
		else:
			raise ValueError('More than one structure in result!')

	@property
	def group_name(self):
		if len(self.group_names) == 1:
			return self.group_names[0]
		else:
			raise ValueError('More than one group name in result!')

	def add(self, *frame_processors):
		new_values = []

		for frame_processor in frame_processors:
			has_structure = frame_processor.structure in self.structures
			has_group_name = frame_processor.group_name in self.group_names

			if not has_structure:
				self.structures.append(frame_processor.structure)

			if not has_group_name:
				self.group_names.append(frame_processor.group_name)

			new_values.extend(frame_processor.values)

		if self.values.ndim > 1:
			self.values = np.vstack([self.values, new_values])
		else:
			self.values = np.array(new_values)


class ProcessorGroup:
	def __init__(self, FrameProcessorObject, structure, atom_groups, *args, **kwargs):
		self.FrameProcessorObject = FrameProcessorObject
		self.structure = structure
		self.atom_groups = self.flatten_groups(atom_groups)
		self.root_dirpath = kwargs.get('root_dirpath', None)

		self.frame_processors = {}
		for group_name, atom_group in self.atom_groups.items():
			if self.root_dirpath is not None:
				filepath = os.path.join(self.root_dirpath, self.get_group_filename(group_name)+'.tsv')
			else:
				filepath = None
			self.frame_processors[group_name] = self.FrameProcessorObject(
				structure, atom_group,
				group_name=group_name, filepath=filepath,
				*args, **kwargs
			)

	def __iter__(self):
		return iter(self.frame_processors.values())

	def get_group_filename(self, group_name):
		if type(group_name) == str:
			return self.structure.label+'_'+group_name
		else:
			return self.structure.label+'_'+'_'.join(group_name)

	def initialize(self):
		for frame_processor in self.frame_processors.values():
			frame_processor.initialize()

	def do(self):
		for frame_processor in self.frame_processors.values():
			frame_processor.do()

	def write(self):
		for frame_processor in self.frame_processors.values():
			frame_processor.write()

	def finish(self):
		for frame_processor in self.frame_processors.values():
			frame_processor.finish()

	@classmethod
	def flatten_groups(cls, atom_groups, parent_key='all'):
		items = {}

		if isinstance(atom_groups, dict):
			for key, value in atom_groups.items():
				if parent_key == 'all':
					new_key = key
				else:
					if isinstance(parent_key, tuple):
						new_key = parent_key + (key,)
					else:
						new_key = (parent_key, key)

				results = cls.flatten_groups(value, new_key)
				items.update(results)
		else:
			items[parent_key] = atom_groups

		return items


class FrameProcessorContainer:
	def __init__(self, FrameProcessorObject, dirpath=None):
		self.FrameProcessorObject = FrameProcessorObject
		self.groups = {}

		self.dirpath = dirpath
		if self.dirpath is not None and not os.path.isdir(self.dirpath):
			os.mkdir(self.dirpath)

	def __iter__(self):
		return iter(groups.values())

	def get_results(self):
		for group in self.groups.values():
			for frame_processor in group:
				yield ProcessorResult(frame_processor)

	def get_result_group_by_class(self):
		results = {}
		for structure, group in self.groups.items():
			for group_name, frame_processor in group.frame_processors.items():
				class_name = frame_processor.structure.class_name

				if (class_name, group_name) in results:
					result = results[(class_name, group_name)]
				else:
					result = ProcessorResult()
					results[class_name, group_name] = result

				result.add(frame_processor)

		return results.values()

	def get_result_by_group_name(self, group_name):
		result = ProcessorResult()
		for structure, group in self.groups.items():
			for current_group_name, frame_processor in group.items():
				if current_group_name == group_name:
					results.add(frame_processor)

		return result

	def get_results_by_group(self):
		results = {}
		for structure, group in self.groups.items():
			for group_name, frame_processor in group.frame_processors.items():
				if group_name in results:
					result = results[group_name]
				else:
					result = ProcessorResult()
					results[group_name] = result

				result.add(frame_processor)

		return results

	def add(self, basic_frame, atom_groups, *args, **kwargs):
		self.groups[basic_frame] = ProcessorGroup(self.FrameProcessorObject, basic_frame, atom_groups, root_dirpath=self.dirpath, *args, **kwargs)

	def initialize(self, basic_frame):
		self.groups[basic_frame].initialize()

	def do(self, basic_frame):
		self.groups[basic_frame].do()

	def write(self, basic_frame):
		self.groups[basic_frame].write()

	def finish(self, basic_frame):
		self.groups[basic_frame].finish()


class ProcedureModule:
	__ID = 0

	def __init__(self, *args, **kwargs):
		self.root_dirpath = kwargs.get('root_dirpath', None)
		if self.root_dirpath is not None:
			self.procedure_dirpath = os.path.join(self.root_dirpath, self._LABEL)
			if not os.path.isdir(self.procedure_dirpath):
				os.mkdir(self.procedure_dirpath)
		else:
			self.procedure_dirpath = None

		self.__id = self.__ID
		self.__ID += 1

		self.dirname = kwargs.get('dirname', '{:0>4d}'.format(self.__id))

	def get_container_dirpath(self, container_dirname):
		if self.procedure_dirpath == None:
			return None
		else:
			return os.path.join(self.procedure_dirpath, container_dirname)

	def execute(self, basic_frame):
		raise NotImplementedError

	def save(self, dirpath):
		raise NotImplementedError

	def terminate(self, basic_frame):
		raise NotImplementedError


class TrajectoryAlignment(ProcedureModule):
	__ID = 0
	_NAME = 'Trajectory alignment'
	_LABEL = 'align'
	_ATOM_SELECTION = 'name CA'

	def __init__(self, structures, reference, *args, **kwargs):
		super().__init__(*args, **kwargs)

		target_atoms = reference.get_residues_by_group('alignment').atoms.select_atoms(self._ATOM_SELECTION)
		self.alignments = FrameProcessorContainer(Alignment, self.get_container_dirpath('rmsds'))
		for structure in structures:
			selected_atoms = structure.get_residues_by_group('alignment').atoms.select_atoms(self._ATOM_SELECTION)
			if len(selected_atoms) == len(target_atoms):
				self.alignments.add(structure, selected_atoms, target_atoms)
			else:
				try:
					logger.error(f'The number of atoms between reference and {structure.label} is different. The superposition may failed and it is due to missing atoms in structure.')
					# Backup plan
					target_residue_indices = reference.domains['alignment'].get_residues().indices
					selected_residue_indices = structure.domains['alignment'].get_residues().indices

					selected_residue_index_iterator = iter(selected_residue_indices)
					selected_residue_index = next(selected_residue_index_iterator)

					keep_selected_atom_indices = []
					keep_target_atom_indices = []

					for target_residue_index in target_residue_indices:
						target_residue = reference.protein.residues[target_residue_index]

						target_atom = target_residue.atoms.select_atoms(self._ATOM_SELECTION)
						if not len(target_atom):
							selected_residue_index = next(selected_residue_index_iterator)
							continue

						selected_residue = structure.protein.residues[selected_residue_index]
						try:
							while selected_residue.resname != target_residue.resname:
								selected_residue_index = next(selected_residue_index_iterator)
								selected_residue = structure.protein.residues[selected_residue_index]

							selected_atom = selected_residue.atoms.select_atoms(self._ATOM_SELECTION)
							if not len(selected_atom):
								continue
							else:
								# Check if previous and next residues are the same
								try:
									is_previous_residue = structure.protein.residues[selected_residue_index-1].resname == reference.protein.residues[target_residue_index-1].resname
									is_next_residue = structure.protein.residues[selected_residue_index+1].resname == reference.protein.residues[target_residue_index+1].resname

									if is_previous_residue and is_next_residue:
										keep_selected_atom_indices.append(selected_atom[0].index)
										keep_target_atom_indices.append(target_atom[0].index)
								except IndexError:
									pass
						except StopIteration:
							break

					keep_selected_atom_indices = np.array(keep_selected_atom_indices)
					keep_target_atom_indices = np.array(keep_target_atom_indices)

					selected_atoms = structure.protein[keep_selected_atom_indices]
					target_atoms = reference.protein[keep_target_atom_indices]

					self.alignments.add(structure, selected_atoms, target_atoms)
				except Exception:
					# Let the program failed to superpose
					logger.error('Backup plan for alignment failed!')
					target_atoms = reference.get_residues_by_group('alignment').atoms.select_atoms(self._ATOM_SELECTION)
					selected_atoms = structure.get_residues_by_group('alignment').atoms.select_atoms(self._ATOM_SELECTION)
					self.alignments.add(structure, selected_atoms, target_atoms)


	def start(self, structure):
		self.alignments.initialize(structure)

	def execute(self, structure):
		self.alignments.do(structure)

	def save(self, structure):
		self.alignments.write(structure)

	def terminate(self, structure):
		self.alignments.finish(structure)


class TrajectoryWriter(ProcedureModule):
	__ID = 0
	_NAME = 'Trajectory writer'
	_LABEL = 'trajwrite'

	def __init__(self, structures, *args, **kwargs):
		super().__init__(*args, **kwargs)

		self.protein_writers = FrameProcessorContainer(FrameWriter, self.procedure_dirpath)
		for structure in structures:
			traj_filepath = os.path.join(self.procedure_dirpath, structure.label)
			self.protein_writers.add(structure, structure.protein, traj_filepath)

	def start(self, structure):
		self.protein_writers.initialize(structure)

	def execute(self, structure):
		self.protein_writers.do(structure)

	def save(self, structure):
		self.protein_writers.write(structure)

	def terminate(self, structure):
		self.protein_writers.finish(structure)    


class CoordinateDeviation(ProcedureModule):
	__ID = 0
	_NAME = 'Coordinate deviation'
	_LABEL = 'deviation'

	def __init__(self, structures):
		super().__init__(*args, **kwargs)

		self.coordinates = FrameProcessorContainer(FrameCoordinate, self.procedure_dirpath)
		for structure in structures:
			self.coordinates.add(structure, structure.protein, dirpath=self.procedure_dirpath)

	def start(self, structure):
		self.coordinates.initialize(structure)

	def execute(self, structure):
		self.coordinates.do(structure)

	def execute(self, structure):
		self.coordinates.save(structure)

	def terminate(self, structure):
		self.coordinates.finish(structure)


class Phelix(ProcedureModule):
	_NAME = 'Phelix'
	_LABEL = 'phelix'
	_ATOM_SELECTION = 'name CA'
	_NEIGHBORS = 2

	def __init__(self, structures, *args, **kwargs):
		super().__init__(*args, **kwargs)

		self.centroids = FrameProcessorContainer(FrameCentroid, self.get_container_dirpath('centroids'))
		for structure in structures:
			transmembrane_extremity_groups = structure.get_transmembrane_extremity_residues(self._NEIGHBORS, self._ATOM_SELECTION)
			self.centroids.add(structure, transmembrane_extremity_groups)

	def start(self, structure):
		self.centroids.initialize(structure)

	def execute(self, structure):
		self.centroids.do(structure)

	def save(self, structure):
		self.centroids.write(structure)

	def terminate(self, structure):
		self.centroids.finish(structure)


class Analyzer:
	_ANALYSIS_OBJECTS = {
		TrajectoryAlignment._LABEL: TrajectoryAlignment,
		TrajectoryWriter._LABEL: TrajectoryWriter,
		CoordinateDeviation._LABEL: CoordinateDeviation,
		Phelix._LABEL: Phelix,
	}

	def __init__(self, structures, reference=None, output_dirpath=None):
		self.reference = reference

		self.structures = []
		for structure in structures:
			if isinstance(structure, frame.Basic):
				self.structures.append(structure)
			else:
				raise TypeError(f'Structure attribute must be {str(frame.BasicFrame)} not {str(type(structure))}')

		self.output_dirpath = output_dirpath
		self.analysis_procedures = []

	def add_analysis(self, *analysis_names, **kwargs):
		for analysis_name in analysis_names:
			AnalysisObject = self._ANALYSIS_OBJECTS[analysis_name]
			self.analysis_procedures.append(AnalysisObject(self.structures, self.reference, root_dirpath=self.output_dirpath, **kwargs))

	def run(self):
		for structure in self.structures:
			for procedure in self.analysis_procedures:
				procedure.start(structure)

			for t, pt in enumerate(structure.trajectory):
				for procedure in self.analysis_procedures:
					procedure.execute(structure)
					procedure.save(structure)

			for procedure in self.analysis_procedures:
				procedure.terminate(structure)