import java.io.PrintWriter;
import org.junit.platform.launcher.Launcher;
import org.junit.platform.launcher.LauncherDiscoveryRequest;
import org.junit.platform.launcher.core.LauncherDiscoveryRequestBuilder;
import org.junit.platform.launcher.core.LauncherFactory;
import org.junit.platform.launcher.listeners.SummaryGeneratingListener;
import static org.junit.platform.engine.discovery.DiscoverySelectors.selectClass;

/** Independent launcher for owned tests only; does not use Maven/target. */
public final class NumericJUnit {
  public static void main(String[] args) {
    LauncherDiscoveryRequest request = LauncherDiscoveryRequestBuilder.request()
        .selectors(selectClass("trimal.util.CppFormatTest"), selectClass("trimal.util.UtilsQuicksortTest"),
                   selectClass("trimal.util.UtilsTest"), selectClass("trimal.stat.SimilarityMatrixTest"),
                   selectClass("trimal.stat.StatisticsGapsTest"), selectClass("trimal.stat.StatisticsConservationTest"))
        .build();
    SummaryGeneratingListener listener = new SummaryGeneratingListener();
    Launcher launcher = LauncherFactory.create(); launcher.registerTestExecutionListeners(listener); launcher.execute(request);
    listener.getSummary().printTo(new PrintWriter(System.out));
    listener.getSummary().printFailuresTo(new PrintWriter(System.out));
    if (listener.getSummary().getTotalFailureCount() != 0) System.exit(1);
  }
}
