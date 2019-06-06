import pyrosetta
import pyrosetta.rosetta
from pyrosetta.rosetta import *

class rotation_mover(pyrosetta.rosetta.protocols.moves.Mover):
	'''A mover that increments along a particular axis'''
	def __init__(self, pose, axis, rotation_mag):
			
		# set active jump info
		self.active_jump = pose.jump(1)
		self.stored_upstream_stub = pose.conformation().upstream_jump_stub(1)
		self.stored_base_centroid = core.chemical.rna.get_rna_base_centroid(pose.residue(2))
		
		# other parameters for mover
		self.axis = axis
		self.rotation_mag = rotation_mag
		
		# generate roation matrix
		self.set_rotation_matrix()	
		
	def get_name(self):
		'''Return name of class.'''
		return self.__class__.__name__
	
	def apply(self, pose):
		
		active_jump = pose.jump(1)
		active_jump.rotation_by_matrix(self.stored_upstream_stub, self.stored_base_centroid, self.rot_matrix)
		pose.set_jump(1, active_jump) 
		
	def set_axis(self, axis):
	
		self.axis = axis
		self.set_rotation_matrix()
	
	def set_rotation_mag(self, rotation_mag):
		
		self.rotation_mag = rotation_mag
		self.set_rotation_matrix()
		
	def set_rotation_matrix(self):
	
		self.rot_matrix = numeric.rotation_matrix(self.axis, self.rotation_mag)	
		
class translation_mover(pyrosetta.rosetta.protocols.moves.Mover):
	'''A mover that increments along a particular axis'''
	def __init__(self, pose, axis, translation_mag):
		
		# set active jump info
		self.active_jump = pose.jump(1)
		
		self.axis = axis
		self.translation_mag = translation_mag
		
		# generate translation vector
		self.set_translation_vector()
		
	def get_name(self):
		'''Return name of class.'''
		return self.__class__.__name__
	
	def apply(self, pose):
		
		active_jump = pose.jump(1)
		active_jump.set_translation( active_jump.get_translation() + self.trans_vector )
		pose.set_jump(1, active_jump)
		
	def set_axis(self, axis):
	
		self.axis = axis
		self.set_translation_vector()
		
	def set_translation_mag(self, translation_mag):

		self.translation_mag = translation_mag
		self.set_translation_vector()
	
	def set_translation_vector(self):

		self.trans_vector = numeric.xyzVector_double_t()
		self.trans_vector.assign(self.axis)
		self.trans_vector.x *= self.translation_mag
		self.trans_vector.y *= self.translation_mag
		self.trans_vector.z *= self.translation_mag