from scipy.integrate import odeint
import pandas as pd

class Reaction:
    def __init__(self, reactions, fixed_concentrations=None) -> None:
        """
        Initialize the Reaction class with a list of reactions.
        
        Parameters:
        reactions (list): A list of reactions with associated rate constants.
        """
        self.reactions = reactions
        self.fixed_concentrations = fixed_concentrations or []
    
    def parse_reaction(self, reaction_str):
        # Парсим реакцию и возвращаем реагенты и продукты
        is_back = False
        if '<->' in reaction_str:
            reactants_str, products_str = reaction_str.split('<->')
            is_back = True
        else:
            reactants_str, products_str = reaction_str.split('->')
        
        reactants = [r.strip() for r in reactants_str.split('+')]
        products = [p.strip() for p in products_str.split('+')]
        
        return reactants, products, is_back

    def create_reaction_dict(self):
        result = {}
        
        for reaction in self.reactions:
            reactants, products, is_back = self.parse_reaction(reaction[0])
            
            # Обработка прямой реакции
            for reactant in reactants:
                if reactant not in result:
                    result[reactant] = []
                result[reactant].append((tuple(reactants), -reaction[1]))
                if is_back:
                    result[reactant].append((tuple(products), reaction[2]))
            
            for product in products:
                if product not in result:
                    result[product] = []
                result[product].append((tuple(reactants), reaction[1]))
                if is_back:
                    result[product].append((tuple(products), -reaction[2]))

        self.reaction_dict = result          


    def rate_of_change(self, concs, t):
        """
        Calculate the rate of change of each substance at a given time point.

        Parameters:
        concs (list): List of current concentrations of the substances.
        t (float): Current time point.

        Returns:
        list: List of rate of changes of the substances.
        """
        conc_dict = {key: value for key, value in zip(self.sorted_keys, concs)}
        rate_dict = {key: 0 for key in self.sorted_keys}


        for substance, reactions_for_substance in self.reaction_dict.items():
            if substance in self.fixed_concentrations:
                rate_dict[substance] = 0
                continue

            for reaction in reactions_for_substance:
                rate = reaction[1]
                reactants = reaction[0]
                for reactant in reactants:
                    if reactant in conc_dict:
                        rate *= conc_dict[reactant]
                    else:
                        raise KeyError(f"Reactant {reactant} not found in the concentration dictionary")
                rate_dict[substance] += rate


        return [rate_dict[key] for key in self.sorted_keys]


    def set_initial_concentrations(self, initial_concs_dict):
        return [initial_concs_dict.get(key, 0) for key in self.sorted_keys]

    def solve(self, initial_concs_dict, time_points):
        # Getting sorted keys from the reaction dictionary
        if not hasattr(self, 'reaction_dict'):
            self.create_reaction_dict()
        if not hasattr(self, 'sorted_keys'):
            self.sorted_keys = sorted(self.reaction_dict.keys())

        # Preparing a list of initial concentrations in the same order as sorted keys
        initial_concs = self.set_initial_concentrations(initial_concs_dict)
        
        # Solving the differential equations
        solution = odeint(self.rate_of_change, initial_concs, time_points)

        # Creating a DataFrame from the solution
        df = pd.DataFrame(solution, columns=self.sorted_keys)
        df.index = time_points

        return df
