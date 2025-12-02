from abc import ABC, abstractmethod

class BasePathwayAPI(ABC):
    @abstractmethod
    def download_pathway_data(self, pathway_id):
        pass

    @abstractmethod
    def download_pathway_image(self, pathway_id):
        pass

    @abstractmethod
    def parse_pathway(self, file_path):
        """Returns: entries, groups"""
        pass