export interface KinematicTestResult {
    name: string;
    passed: boolean;
    error?: string;
    details?: any;
}
/**
 * Test shape oscillation natural frequencies
 *
 * Verify that shape modes oscillate at their natural frequencies
 */
export declare function testShapeOscillationFrequencies(): KinematicTestResult;
/**
 * Test bubble translation with Bjerknes force
 *
 * Verify that bubble moves in response to acoustic gradient
 */
export declare function testBubbleTranslation(): KinematicTestResult;
/**
 * Test shape-radial coupling
 *
 * Verify that shape oscillations affect radial dynamics
 */
export declare function testShapeRadialCoupling(): KinematicTestResult;
/**
 * Run all kinematic tests
 */
export declare function runAllKinematicTests(): {
    passed: number;
    failed: number;
    results: KinematicTestResult[];
};
//# sourceMappingURL=kinematicTests.d.ts.map