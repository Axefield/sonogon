export interface TestResult {
    name: string;
    passed: boolean;
    error?: string;
    details?: any;
}
/**
 * Test adiabatic compression
 *
 * For adiabatic compression: P * V^gamma = constant
 * Verify that P * V^gamma remains constant during compression
 */
export declare function testAdiabaticCompression(): TestResult;
/**
 * Test Saha equilibrium
 *
 * Verify that electron density matches Saha equation prediction
 * at equilibrium conditions
 */
export declare function testSahaEquilibrium(): TestResult;
/**
 * Test energy conservation
 *
 * Verify that total energy is approximately conserved over a cycle
 */
export declare function testEnergyConservation(): TestResult;
/**
 * Test plasma frequency calculation
 *
 * Verify plasma frequency formula: ωp = sqrt(ne * e² / (ε0 * me))
 */
export declare function testPlasmaFrequency(): TestResult;
/**
 * Run all validation tests
 */
export declare function runAllTests(): {
    passed: number;
    failed: number;
    results: TestResult[];
};
export { testAdiabaticScaling, testMinnaertFrequency, testEnergyBudgetClosure, testPlasmaEquilibrium, runAllPhysicsValidationTests, } from "./physicsValidation";
//# sourceMappingURL=tests.d.ts.map