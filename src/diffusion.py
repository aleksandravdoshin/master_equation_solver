from reactor import Reaction
import duckdb
connection = duckdb.connect(database=':memory:', read_only=False)

class Reaction_Diffusion(Reaction):
    def __init__(self, reactions, capacity=None, fixed_concentrations=None) -> None:
        super().__init__(reactions, fixed_concentrations)
        self.capacity = capacity
        if capacity is not None:
            self._normalize_capacity()

    def parse_diffusion(self, reaction_str):
        # Парсим реакцию и возвращаем реагенты и продукты
        is_back = False
        if '<D>' in reaction_str:
            reactants_str, products_str = reaction_str.split('<D>')
            is_back = True
        else:
            reactants_str, products_str = reaction_str.split('D>')
        
        reactants = [r.strip() for r in reactants_str.split('+')]
        products = [p.strip() for p in products_str.split('+')]
        
        return reactants, products, is_back
    
    def _normalize_capacity(self):
        norm_factor = sum(self.capacity.values())
        for key in self.capacity:
            self.capacity[key] = self.capacity[key] / norm_factor

    def parse_reaction(self, reaction_str):
        if ('<->' in reaction_str) or ('->' in reaction_str):
            return super().parse_reaction(reaction_str), False
        else:
            return self.parse_diffusion(reaction_str), True

    def create_reaction_dict(self):
        result = {}
        
        for reaction in self.reactions:
            (reactants, products, is_back), is_diffusion = self.parse_reaction(reaction[0])
            # Обработка прямой реакции
            for reactant in reactants:
                if reactant not in result:
                    result[reactant] = []
                coefficient = 1/self.capacity[reactant] if is_diffusion else 1
                result[reactant].append((tuple(reactants), -reaction[1] * coefficient))
                if is_back:
                    coefficient = 1/self.capacity[products[0]] if is_diffusion else 1
                    result[reactant].append((tuple(products), reaction[2] * coefficient))
            for product in products:
                if product not in result:
                    result[product] = []
                
                coefficient = 1/self.capacity[reactant] if is_diffusion else 1
                result[product].append((tuple(reactants), reaction[1] * coefficient))
                if is_back:
                    coefficient = 1/self.capacity[product] if is_diffusion else 1
                    result[product].append((tuple(products), -reaction[2] * coefficient))

        self.reaction_dict = result  
